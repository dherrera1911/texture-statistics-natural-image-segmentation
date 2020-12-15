library(tidyr)
library(dplyr)
library(glmnet)
library(caret)
library(vcd)
source("./analysis_functions.R")
set.seed(2691)

segmentStatsFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
agreementFile <- "../../data/BSD_results/6_classification_agreement.Rds"
experimentTexturesFile <- "../../data/experiment_stats/experiment_textures_stats.csv"
experimentResultsFile <- "../../data/experiment_stats/experimentResults.RDS"
saveBetterPredictions <- "../../data/BSD_results/9_better_stats_predictions.RDS"
saveSegmentationPredictions <- "../../data/BSD_results/9_segmentation_predictions.RDS"

repRuns <- 20

#############################
# load data
#############################
segmentStats <- read.csv(segmentStatsFile) %>%
  make_task_BSD(.)
agreementDf <- readRDS(agreementFile) %>%
  as_tibble(.) %>%
  remove_constant_stats(.)
expTextStats <- read.csv(experimentTexturesFile) %>%
  remove_constant_stats(.) %>%
  as_tibble(.)
expResults <- readRDS(experimentResultsFile)

################################################
# Test whether the segmentation model classifies the
# experiment textures as belonging to the same or to
# different segments
################################################

statisticsNames <- get_statistics_names(names(expTextStats))
designNamesExp <- statisticsNames[["design"]]
statisticsNames <- statisticsNames[which(names(statisticsNames)!="design")]

statsTypes <- list("FAS", "HOS", c("FAS", "HOS"))
resultsDf <- data.frame(texture=character(), rep=integer(),
                        FAS=integer(), HOS=integer(),
                        FASHOS=integer())

expPairsStats <- group_by(expTextStats, texture) %>%
  dplyr::mutate_each(funs(abs(. - .[type=="original"])), -type) %>%
  dplyr::filter(., type!="original") %>%
  dplyr::select(., -type) %>%
  dplyr::mutate(., same=0) %>%
  ungroup(.)

for (l in c(1:repRuns)) {
  tempDf <- data.frame(texture=as.character(expPairsStats$texture),
                       rep=l, same_FAS=0, same_HOS=0, same_FASHOS=0,
                       activation_FAS=0, activation_HOS=0,
                       activation_FASHOS=0)
  for (m in c(1:length(statsTypes))) {
    statsToUse <- get_statistics_names(names(expPairsStats))[statsTypes[[m]]] %>%
      unlist_names(.)
    expSegment <- train_test_ridge(trainData=segmentStats,
                                         testData=expPairsStats,
                                         statsToUse=statsToUse,
                                         balanceWeights=TRUE,
                                         subsetsPCA=NA,
                                         labelColumn="same")
    statsTypeStr <- paste(statsTypes[[m]], collapse="")
    activationCol <- paste("activation_", statsTypeStr, sep="")
    predictionCol <- paste("same_", statsTypeStr, sep="")
    tempDf[[activationCol]] <- as.numeric(expSegment$modelResponses)
    tempDf[[predictionCol]] <- as.integer(expSegment$predictions)
  }
  print(paste("Progress:  ", l, "/", repRuns, sep=""))
  resultsDf <- rbind(tempDf)
  saveRDS(resultsDf, saveSegmentationPredictions)
}

################################################
# Train the model on all the data and test on the textures
# from the original experiment
################################################

expPairsStats$betterHOS <- factor(rep("0",4), levels=c("0", "1"))
expPairsStats$betterFASHOS <- factor(rep("0",4), levels=c("0", "1"))

statsToUse <- get_statistics_names(names(expTextStats))[c("FAS", "HOS")] %>%
  unlist_names(.)

predictBetterStats <- data.frame(texture=character(),
                                 betterHOSActivation=numeric(),
                                 betterHOS=integer(),
                                 betterFASHOSActivation=numeric(),
                                 betterFASHOS=integer(),
                                 rep=integer())

for (l in c(1:repRuns)) {
  expPredHOS <- train_test_ridge(trainData=agreementDf,
                                       testData=expPairsStats,
                                       statsToUse=statsToUse,
                                       balanceWeights=TRUE,
                                       subsetsPCA=NA,
                                       labelColumn="betterHOS")
  expPredFASHOS <- train_test_ridge(trainData=agreementDf, 
                                          testData=expPairsStats,
                                          statsToUse=statsToUse,
                                          balanceWeights=TRUE,
                                          subsetsPCA=NA,
                                          labelColumn="betterFASHOS")
  tempDf <- data.frame(texture=expPairsStats$texture,
                       betterHOSActivation=as.numeric(expPredHOS$modelResponses),
                       betterHOS=as.integer(expPredHOS$predictions),
                       betterFASHOSActivation=as.numeric(expPredFASHOS$modelResponses),
                       betterFASHOS=as.integer(expPredFASHOS$predictions),
                       rep=l)
  predictBetterStats <- rbind(predictBetterStats, tempDf)
  saveRDS(predictBetterStats, saveBetterPredictions)
  print(paste("Progress:  ", l, "/", repRuns, sep=""))
}


