library(tidyr)
library(dplyr)
library(glmnet)
source("./analysis_functions.R")
set.seed(2691)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_pixNorm.csv"
dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"

saveResults <- "../../data/BSD_results/BSD_results.Rds"
saveAgreementResults <- "../../data/BSD_results/BSD_agreement.Rds"
saveCoefs <- "../../data/BSD_results/BSD_coefs.Rds"

alpha <- 0
usePixelInfo <- TRUE
repExp <- 20

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

# generate lists and dfs
resultsDf <- data.frame(usePix = logical(), Pix = double(),
                        FA = double(), HOS = double(), FA_HOS = double())
agreementDf <- data.frame(usePix = logical(), FAS_HOS = double(),
                        FAS_FASHOS = double())
modelFitPix <- list()
modelFitFAS <- list()
modelFitHOS <- list()
modelFitFASHOS <- list()

allDataTask <- make_task_BSD(segmentStats)

for (usePixelInfo in c(TRUE, FALSE)) {
  for (r in 1:repExp) {
    nSegments <- length(unique(segmentStats$ImageName))
    sampleSegments <- sample(unique(segmentStats$ImageName))
    trainSegments <- sampleSegments[1:floor(nSegments*4/5)]
    testSegments <- sampleSegments[ceiling(nSegments*4/5):nSegments]

    trainData <- dplyr::filter(allDataTask, ImageName %in% trainSegments) %>%
      droplevels(.)
    testData <- dplyr::filter(allDataTask, ImageName %in% testSegments) %>%
      droplevels(.)

    # calculate weights to even out classes
    ratioSame <- sum(trainData$same) / sum(1-trainData$same)
    weights <- rep(1, nrow(trainData))
    weights[trainData$same == 1] <- (1/ratioSame)

    # extract labels
    trainLabel <- trainData$same
    testLabel <- testData$same

    # model with only Pixels
    if (usePixelInfo) {
      trainStatsPix <- dplyr::select(trainData, all_of(namesPix))
      testStatsPix <- dplyr::select(testData, all_of(namesPix))
      modelFitPix[[r]] <- cv.glmnet(x = as.matrix(trainStatsPix),
                           y = trainLabel,
                           weights = weights,
                           family="binomial",
                           alpha = alpha,
                           standardize=TRUE,
                           type.measure="class")
      modelPredictionsPix <- predict(modelFitPix[[r]], as.matrix(testStatsPix),
                                     type = "class")
      predictionOutcomePix <- mean(as.integer(modelPredictionsPix == testLabel))
    } else {
      modelFitPix[[r]] <- NA
      predictionOutcomePix <- NA
    }

    # extract the statistics to use in each model
    if (usePixelInfo) {
      trainStatsFAS <- dplyr::select(trainData, all_of(c(namesPix, namesFAS)))
      testStatsFAS <- dplyr::select(testData, all_of(c(namesPix, namesFAS)))
      trainStatsHOS <- dplyr::select(trainData, all_of(c(namesPix, namesHOS)))
      testStatsHOS <- dplyr::select(testData, all_of(c(namesPix, namesHOS)))
      trainStatsFASHOS <- dplyr::select(trainData, all_of(c(namesPix, namesFAS, namesHOS)))
      testStatsFASHOS <- dplyr::select(testData, all_of(c(namesPix, namesFAS, namesHOS)))
    } else {
      trainStatsFAS <- dplyr::select(trainData, all_of(namesFAS))
      testStatsFAS <- dplyr::select(testData, all_of(namesFAS))
      trainStatsHOS <- dplyr::select(trainData, all_of(namesHOS))
      testStatsHOS <- dplyr::select(testData, all_of(namesHOS))
      trainStatsFASHOS <- dplyr::select(trainData, all_of(c(namesFAS, namesHOS)))
      testStatsFASHOS <- dplyr::select(testData, all_of(c(namesFAS, namesHOS)))
    }

    # FAS stats model
    modelFitFAS[[r]] <- cv.glmnet(x = as.matrix(trainStatsFAS),
                         y = trainLabel,
                         weights = weights,
                         family="binomial",
                         alpha = alpha,
                         standardize=TRUE,
                         type.measure="class")
    modelPredictionsFAS <- predict(modelFitFAS[[r]], as.matrix(testStatsFAS),
                                   type = "class")
    predictionOutcomeFAS <- mean(as.integer(modelPredictionsFAS == testLabel))

    # HOS stats model
    modelFitHOS[[r]] <- cv.glmnet(x = as.matrix(trainStatsHOS),
                         y = trainLabel,
                         weights = weights,
                         family="binomial",
                         alpha = alpha,
                         standardize=TRUE,
                         type.measure="class")
    modelPredictionsHOS <- predict(modelFitHOS[[r]], as.matrix(testStatsHOS),
                                   type = "class")
    predictionOutcomeHOS <- mean(as.integer(modelPredictionsHOS == testLabel))

    # all stats model
    modelFitFASHOS[[r]] <- cv.glmnet(x = as.matrix(trainStatsFASHOS),
                         y = trainLabel,
                         weights = weights,
                         family="binomial",
                         alpha = alpha,
                         standardize=TRUE,
                         type.measure="class")
    modelPredictionsFASHOS <- predict(modelFitFASHOS[[r]], as.matrix(testStatsFASHOS),
                                      type = "class")
    predictionOutcomeFASHOS <- mean(as.integer(modelPredictionsFASHOS == testLabel))
    
    results <- c(usePixelInfo, predictionOutcomePix, predictionOutcomeFAS,
                 predictionOutcomeHOS, predictionOutcomeFASHOS)
    resultsDf[nrow(resultsDf)+1,] <- results
  }
}

# save results
saveRDS(resultsDf, saveResults)
saveRDS(agreementDf, saveAgreementResults)

