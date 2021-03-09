library(tidyr)
library(dplyr)
library(glmnet)
source("./analysis_functions.R")
set.seed(2691)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_pixNorm.csv"
dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
savePerformanceResults <- "../../data/BSD_results/8_subset_stats_performance.RDS"
saveTexturePredictions <- "../../data/BSD_results/8_texture_predictions_subsets.RDS"
dataSplits <- 10

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
statisticsNames <- get_statistics_names(parNames, subsetHOS=TRUE)
designNames <- statisticsNames$design
statisticsNames <- statisticsNames[which(names(statisticsNames)!="design")]

#############################
#generate template of design matrix for one repetition
#############################
HOStypes <- names(statisticsNames$HOS)
FAS <- c(0)
acm <- c(0,1)
cmc <- c(0,1)
pmc <- c(0,1)
prc <- c(0,1)

designMatrixTemp <- expand.grid(FAS, acm, cmc, pmc, prc) %>%
  dplyr::mutate(., rep=NA, performance=NA) %>%
  dplyr::rename(., FAS=Var1, acm=Var2, cmc=Var3,
                pmc=Var4, prc=Var5) %>%
  dplyr::filter(., (FAS+acm+cmc+pmc+prc)!=0)

modelPerformance <- NULL
texturePredictions <- NULL


#################################################
# Make test splits to cover all images
#################################################
nSegments <- length(unique(segmentStats$ImageName))
sampleSegments <- sample(unique(segmentStats$ImageName))
splitSize <- ceiling(nSegments/dataSplits)
testSegmentsList <- split(sampleSegments,
                          ceiling(seq_along(sampleSegments)/splitSize))

#############################
#Put together the pairs of patches
#############################
allDataTask <- make_task_BSD(segmentStats)

#############################
#Fit the models
#############################
for (r in 1:dataSplits) {
  # get the segments to test, and split the rest into two training sets
  testSegments <- testSegmentsList[[r]]
  trainSegments <- sampleSegments[which(!sampleSegments %in% testSegments)]

  trainData <- dplyr::filter(allDataTask, ImageName %in% trainSegments) %>%
    droplevels(.)
  testData <- dplyr::filter(allDataTask, ImageName %in% testSegments) %>%
    droplevels(.)
  
  copyTemplate <- designMatrixTemp %>%
    dplyr::mutate(., rep = r)

  textureResultsDf <- dplyr::select(testData, same, designNames)

  print(paste("Split: ", r, "/", dataSplits, sep=""))

  for (m in c(1:nrow(copyTemplate))) {
    statsInd <- which(c(copyTemplate[m,c("acm", "cmc", "pmc", "prc")])==1)
    trialHOS <- HOStypes[statsInd]
    if (length(trialHOS)>0) {
      trialHOSList <- statisticsNames$HOS[trialHOS]
      trialStatsVec <- unlist_names(trialHOSList) 
      modelName <- paste(trialHOS, collapse="_")
    } else {
      trialStatsVec <- NULL
      modelName <- NULL
    }
    if (copyTemplate$FAS[m] == 1) {
      trialStatsVec <- c(statisticsNames$FAS, trialStatsVec)
      if (length(trialHOS) > 0) {
        modelName <- paste("FAS_", modelName, sep="")
      } else {
        modelName <- "FAS"
      }
    }
    print(paste("Row: ", m, "/", nrow(copyTemplate), sep=""))
    modelOutcome <- train_test_ridge(trainData=trainData, testData=testData,
                     statsToUse=trialStatsVec, balanceWeights=TRUE, subsetsPCA=NA)
    copyTemplate$performance[m] <- modelOutcome$accuracy
    textureResultsDf[[modelName]] <- modelOutcome$predictions
  }
  modelPerformance <- rbind(modelPerformance, copyTemplate)
  texturePredictions <- rbind(texturePredictions, textureResultsDf)
  saveRDS(modelPerformance, savePerformanceResults)
  saveRDS(texturePredictions, saveTexturePredictions)
}

