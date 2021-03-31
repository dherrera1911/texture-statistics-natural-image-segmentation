library(tidyr)
library(dplyr)
library(glmnet)
source("./analysis_functions.R")
set.seed(2691)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_pixNorm.csv"
dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
saveResults <- "../../data/BSD_results/3_BSD_results_noNorm.RDS"
repExp <- 20
subsetPCA <- NA
pcaVar <- NA
normalizeData <- FALSE

#############################
# load data
#############################
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
# make name list for doing PCA within subsets of HOS if needed
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
  dplyr::mutate(., rep = NA, performance = NA) %>%
  dplyr::rename(., pixel = Var1, FAS = Var2, HOS = Var3) %>%
  dplyr::filter(., !(pixel==0 & FAS==0 & HOS==0))
resultsDf <- NULL

#############################
#Put together the pairs of patches
#############################
allDataTask <- make_task_BSD(segmentStats)

#############################
#Fit the models
#############################
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
  
  copyTemplate <- designMatrixTemp %>%
    dplyr::mutate(., rep = r)

  for (m in c(1:nrow(copyTemplate))) {
    statsInd <- which(c(copyTemplate[m,c("pixel", "FAS", "HOS")])==1)
    trialTypes <- statsTypes[statsInd]
    trialStatsList <- statisticsNames[trialTypes]
    trialStatsVec <- unlist_names(trialStatsList) 
    if (is.na(subsetPCA)) {
      pcaStatsSubsets <- NA
    } else if (subsetPCA) {
      pcaStatsSubsets <- statisticsNamesFiner
    } else {
      pcaStatsSubsets <- list(all=trialStatsVec)
    }
    modelOutcome <- train_test_ridge(trainData=trainData,
                                     testData=testData,
                                     statsToUse=trialStatsVec,
                                     balanceWeights=TRUE,
                                     subsetsPCA=pcaStatsSubsets,
                                     varianceRetained=pcaVar,
                                     normalizeData=normalizeData)
    copyTemplate$performance[m] <- modelOutcome$accuracy
    print(paste("Rep: ", r,"/", repExp, "     Row: ", m, "/",
                nrow(copyTemplate), sep=""))
  }
  resultsDf <- rbind(resultsDf, copyTemplate)
  saveRDS(resultsDf, saveResults)
}

