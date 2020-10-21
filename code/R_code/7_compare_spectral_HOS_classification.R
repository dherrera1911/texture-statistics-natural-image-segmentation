library(tidyr)
library(dplyr)
library(glmnet)
library(caret)
library(vcd)
source("./analysis_functions.R")
set.seed(2691)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_pixNorm.csv"
dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"

saveBetterHOSpairs <- "../../data/BSD_stats/pairs_HOS_better.Rds"
saveBetterFASHOSpairs <- "../../data/BSD_stats/pairs_FASHOS_better.Rds"

alpha <- 0
dataSplits <- 10

# load data
segmentStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.)

# get the indices of the different types of stats to use
parNames <- names(segmentStats)
statisticsNames <- get_statistics_names(parNames)
namesPix <- statisticsNames$pixel
namesFAS <- statisticsNames$FAS
namesHOS <- statisticsNames$HOS
designNames <- statisticsNames$design

# initialize lists and dfs
resultsDf <- data.frame(FA = double(), HOS = double(), FA_HOS = double())
modelFitFAS <- list()
modelFitHOS <- list()
modelFitFASHOS <- list()
betterHOSDf <- NULL
betterFASHOSDf <- NULL

# make splits of images to test for better HOS
nSegments <- length(unique(segmentStats$ImageName))
sampleSegments <- sample(unique(segmentStats$ImageName))
splitSize <- ceiling(nSegments/dataSplits)
testSegmentsList <- split(sampleSegments,
                          ceiling(seq_along(sampleSegments)/splitSize))

allDataTask <- make_task_BSD(segmentStats)

for (r in 1:dataSplits) {
  # get the segments to test, and split the rest into two training sets
  testSegments <- testSegmentsList[[r]]
  trainSegmentsAll <- sampleSegments[which(!sampleSegments %in% testSegments)]
  trainSegmentsAll <- sample(trainSegmentsAll)
  nTrainSegs <- length(trainSegmentsAll)
  halfTrainSegs <- round(nTrainSegs/2)
  trainSegments1 <- trainSegmentsAll[1:halfTrainSegs]
  trainSegments2 <- trainSegmentsAll[(halfTrainSegs+1):nTrainSegs]

  trainDataSplit <- list()
  trainDataSplit[[1]] <- dplyr::filter(allDataTask, ImageName %in% trainSegments1) %>%
    droplevels(.) %>%
    as_tibble(.)
  trainDataSplit[[2]] <- dplyr::filter(allDataTask, ImageName %in% trainSegments2) %>%
    droplevels(.) %>%
    as_tibble(.)
  testData <- dplyr::filter(allDataTask, ImageName %in% testSegments) %>%
    droplevels(.) %>%
    as_tibble(.)

  # train and test models for the two training set splits
  modelFitFAS[[r]] <- list()
  modelFitHOS[[r]] <- list()
  modelFitFASHOS[[r]] <- list()
  betterHOSPairs <- list()
  betterFASHOSPairs <- list()
  for (s in c(1:2)) {
    trainData <- trainDataSplit[[s]]
    # calculate weights to even out classes
    ratioSame <- sum(trainDataSplit[[s]]$same) /
      sum(1-trainDataSplit[[s]]$same)
    weights <- rep(1, nrow(trainData))
    weights[trainData$same == 1] <- (1/ratioSame)

    # extract labels
    trainLabel <- trainData$same
    testLabel <- testData$same

    trainStatsFAS <- dplyr::select(trainData, all_of(namesFAS))
    testStatsFAS <- dplyr::select(testData, all_of(namesFAS))
    trainStatsHOS <- dplyr::select(trainData, all_of(namesHOS))
    testStatsHOS <- dplyr::select(testData, all_of(namesHOS))
    trainStatsFASHOS <- dplyr::select(trainData, all_of(c(namesFAS, namesHOS)))
    testStatsFASHOS <- dplyr::select(testData, all_of(c(namesFAS, namesHOS)))

    # FAS stats model
    modelFitFAS[[r]][[s]] <- cv.glmnet(x = as.matrix(trainStatsFAS),
                         y = trainLabel,
                         weights = weights,
                         family="binomial",
                         alpha = alpha,
                         standardize=TRUE,
                         type.measure="class")
    modelPredictionsFAS <- predict(modelFitFAS[[r]][[s]],
                                   as.matrix(testStatsFAS), type = "class")
    predictionOutcomeFAS <- mean(as.integer(modelPredictionsFAS == testLabel))

    # HOS stats model
    modelFitHOS[[r]][[s]] <- cv.glmnet(x = as.matrix(trainStatsHOS),
                         y = trainLabel,
                         weights = weights,
                         family="binomial",
                         alpha = alpha,
                         standardize=TRUE,
                         type.measure="class")
    modelPredictionsHOS <- predict(modelFitHOS[[r]][[s]],
                                   as.matrix(testStatsHOS), type = "class")
    predictionOutcomeHOS <- mean(as.integer(modelPredictionsHOS == testLabel))

    # all stats model
    modelFitFASHOS[[r]][[s]] <- cv.glmnet(x = as.matrix(trainStatsFASHOS),
                         y = trainLabel,
                         weights = weights,
                         family="binomial",
                         alpha = alpha,
                         standardize=TRUE,
                         type.measure="class")
    modelPredictionsFASHOS <- predict(modelFitFASHOS[[r]][[s]],
                                      as.matrix(testStatsFASHOS), type = "class")
    predictionOutcomeFASHOS <- mean(as.integer(modelPredictionsFASHOS == testLabel))

    ##### Get pairs where FAS failed but HOS or FAS+HOS succeeded
    wrongPairsFAS <- which(as.numeric(modelPredictionsFAS) != testLabel)
    wrongPairsHOS <- which(as.numeric(modelPredictionsHOS) != testLabel)
    wrongPairsFASHOS <- which(as.numeric(modelPredictionsFASHOS) != testLabel)
    correctPairsFAS <- which(as.numeric(modelPredictionsFAS) == testLabel)
    correctPairsHOS <- which(as.numeric(modelPredictionsHOS) == testLabel)
    correctPairsFASHOS <- which(as.numeric(modelPredictionsFASHOS) == testLabel)

    betterHOSPairs[[s]] <- wrongPairsFAS[which(wrongPairsFAS %in%
                                               correctPairsHOS)]
    betterFASHOSPairs[[s]] <- wrongPairsFAS[which(wrongPairsFAS %in%
                                                  correctPairsFASHOS)]
  }

  betterHOSAgree <- betterHOSPairs[[1]][which(betterHOSPairs[[1]] %in%
                                                 betterHOSPairs[[2]])]
  betterFASHOSAgree <- betterFASHOSPairs[[1]][which(betterFASHOSPairs[[1]] %in%
                                                       betterFASHOSPairs[[2]])]
  betterHOSDf <- rbind(betterHOSDf, testData[betterHOSAgree,])
  betterFASHOSDf <- rbind(betterFASHOSDf, testData[betterFASHOSAgree,])
}


# save dataframe with segment pairs where HOS is better
saveRDS(betterHOSDf, saveBetterHOSpairs)
saveRDS(betterFASHOSDf, saveBetterFASHOSpairs)

