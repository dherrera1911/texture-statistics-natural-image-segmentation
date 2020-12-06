library(tidyr)
library(dplyr)
library(glmnet)
source("./analysis_functions.R")
set.seed(2691)

dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
saveResults <- "../../data/BSD_results/4_BSD_dnn.Rds"

repExp <- 10
layerUnits <- list(c(30), c(30, 10), c(30, 10, 2), c(50, 20), c(50, 20, 5), c(10, 2))
regularizationWeight <- c(0.001, 0.002, 0.004)
epochs <- c(150, 200, 300)

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

for (r in 1:repExp) {
  nSegments <- length(unique(segmentStats$ImageName))
  sampleSegments <- sample(unique(segmentStats$ImageName))
  trainSegments <- sampleSegments[1:floor(nSegments*4/5)]
  testSegments <- sampleSegments[(floor(nSegments*4/5)+1):nSegments]

  trainData <- dplyr::filter(allDataTask, ImageName %in% trainSegments) %>%
    droplevels(.)
  testData <- dplyr::filter(allDataTask, ImageName %in% testSegments) %>%
    droplevels(.)

  for (arq in c(1:length(layerUnits))) {
    archString <- paste(layerUnits[[arq]], sep="-", collapse="-")
    for (regW in regularizationWeight) {
      for (ep in epochs) {
        # split data into train and test set
        copyTemplate <- designMatrixTemp %>%
          dplyr::mutate(., rep = r)
        copyTemplate$architecture <- archString
        copyTemplate$regularization <- regW
        copyTemplate$epochs <- ep
        for (m in c(1:nrow(copyTemplate))) {
          statsInd <- which(c(copyTemplate[m,c("pixel", "FAS", "HOS")])==1)
          trialTypes <- statsTypes[statsInd]
          trialStatsList <- statisticsNames[trialTypes]
          trialStatsVec <- unlist_names(trialStatsList) 
          modelOutcome <- train_test_dnn(trainData=trainData, testData=testData,
                           statsToUse=trialStatsVec, balanceWeights=TRUE,
                           subsetsPCA=list(all=trialStatsVec),
                           layerUnits=layerUnits[[arq]],
                           regularizationWeight=regW,
                           epochs=ep)
          copyTemplate$performance[m] <- modelOutcome$accuracy
          progressText <- paste("Rep:", r, "  Architecture:", archString,
                                "  RegW:", regW, "  Epochs:", ep,
                                "  Stats comb:", m, sep="")
          print(progressText)
        }
      resultsDf <- rbind(resultsDf, copyTemplate)
      saveRDS(resultsDf, saveResults)
      }
    }
  }
}

# save results

