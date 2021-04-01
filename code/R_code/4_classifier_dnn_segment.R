library(tidyr)
library(dplyr)
library(glmnet)
source("./analysis_functions.R")
set.seed(2691)

dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
saveResults <- "../../data/BSD_results/4_BSD_dnn_noNorm.Rds"
saveResultsHistory <- "../../data/BSD_results/4_BSD_dnn_history_noNorm.Rds"

repExp <- 20
layerUnits <- c(30)
regularizationWeight <- 0.01
epochs <- 300
subsetPCA <- FALSE #whether to do the PCA separately for each group of statistics
finerSubset <- FALSE # whether to subset the HOS subsets separately
normalizeData <- FALSE

# load data
segmentStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.) %>%
  remove_constant_stats(.)

#############################
# get the names of the different stats to use
#############################
parNames <- names(segmentStats)
statisticsNames <- get_statistics_names(parNames)
designNames <- statisticsNames$design
statisticsNames <- statisticsNames[which(names(statisticsNames)!="design")]
# make name list for doing PCA within subsets of HOS
statisticsNamesFiner <- get_statistics_names(parNames, subsetHOS=TRUE)
statisticsNamesFiner <- statisticsNamesFiner[which(names(statisticsNamesFiner)!="design")]
statisticsNamesFiner <- c(list(pixel=statisticsNamesFiner$pixel),
                          list(FAS=statisticsNamesFiner$FAS),
                          statisticsNamesFiner$HOS)

#############################
#generate template of design matrix for one repetition
#############################
pixel <- c(0,1)
FAS <- c(0,1)
HOS <- c(0,1)
statsTypes <- c("pixel", "FAS", "HOS")
designMatrixTemp <- expand.grid(pixel, FAS, HOS) %>%
  dplyr::mutate(., rep=NA, performance=NA) %>%
  dplyr::rename(., pixel=Var1, FAS=Var2, HOS=Var3) %>%
  dplyr::filter(., !(pixel==0 & FAS==0 & HOS==0))
resultsDf <- NULL

#############################
#Put together the pairs of patches
#############################
allDataTask <- make_task_BSD(segmentStats)
trainingHistory <- list()

for (r in 1:repExp) {
  # split data into train and test set
  nSegments <- length(unique(segmentStats$ImageName))
  sampleSegments <- sample(unique(segmentStats$ImageName))
  trainSegments <- sampleSegments[1:floor(nSegments*4/5)]
  testSegments <- sampleSegments[(floor(nSegments*4/5)+1):nSegments]

  trainData <- dplyr::filter(allDataTask, ImageName %in% trainSegments) %>%
    droplevels(.)
  testData <- dplyr::filter(allDataTask, ImageName %in% testSegments) %>%
    droplevels(.)

  archString <- paste(layerUnits, sep="-", collapse="-")
  copyTemplate <- designMatrixTemp %>%
    dplyr::mutate(., rep=r)
  copyTemplate$architecture <- archString
  copyTemplate$regularization <- regularizationWeight
  copyTemplate$epochs <- epochs
  copyTemplate$PCA_subsetting <- subsetPCA
  for (m in c(1:nrow(copyTemplate))) {
    statsInd <- which(c(copyTemplate[m,c("pixel", "FAS", "HOS")])==1)
    trialTypes <- statsTypes[statsInd]
    trialStatsList <- statisticsNames[trialTypes]
    trialStatsVec <- unlist_names(trialStatsList) 
    if (is.na(subsetPCA)) {
      pcaStatsSubsets <- NA
    } else if (subsetPCA) {
      if (finerSubset) {
        pcaStatsSubsets <- statisticsNamesFiner
      } else {
        pcaStatsSubsets <- trialStatsList
      }
    } else {
      pcaStatsSubsets <- list(all=trialStatsVec)
    }
    modelOutcome <- train_test_dnn(trainData=trainData, testData=testData,
                     statsToUse=trialStatsVec, balanceWeights=TRUE,
                     subsetsPCA=pcaStatsSubsets,
                     layerUnits=layerUnits,
                     regularizationWeight=regularizationWeight,
                     epochs=epochs,
                     normalizeData=normalizeData)
    copyTemplate$performance[m] <- modelOutcome$accuracy
    progressText <- paste("Rep:", r, "  Architecture:", archString,
                          "  RegW:", regularizationWeight, "  Epochs:", epochs,
                          "  Stats comb:", m, sep="")
    # save the training history
    statsName <- assign_stat_name(pixel=copyTemplate[m,"pixel"],
                                 FAS=copyTemplate[m,"FAS"],
                                 HOS=copyTemplate[m,"HOS"])
    trainingHistory[[statsName]] <- rbind(trainingHistory[[statsName]],
                                          modelOutcome$accuracyHistory)
  }
  resultsDf <- rbind(resultsDf, copyTemplate)
  saveRDS(resultsDf, saveResults)
  saveRDS(trainingHistory, saveResultsHistory)
}

