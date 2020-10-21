library(tidyr)
library(dplyr)
library(glmnet)
source("./analysis_functions.R")
set.seed(2691)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_pixNorm.csv"
dataFile <- "../../data/texture_stats/texture_stats_statsNorm.csv"

saveResults <- "../../data/texture_results/texture_results.Rds"
saveAgreementResults <- "../../data/texture_results/texture_agreement.Rds"
saveCoefs <- "../../data/texture_results/texture_coefs.Rds"

nRep <- 5
repExp <- 20
alpha <- 0
usePixelInfo <- TRUE

# load data
textureStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.)

# get the indices of the different types of stats to use
parNames <- names(textureStats)
statisticsNames <- get_statistics_names(parNames)
namesPix <- statisticsNames$pixel
namesFAS <- statisticsNames$FAS
namesHOS <- statisticsNames$HOS
designNames <- statisticsNames$design


# generate lists and dfs
resultsDf <- data.frame(usePix = logical(), Pix = double(),
                        FA = double(), HOS = double(),
                        FA_HOS = double())
agreementDf <- data.frame(usePix = logical(), FAS_HOS = double(),
                        FAS_FASHOS = double())
modelFitPix <- list()
modelFitFAS <- list()
modelFitHOS <- list()
modelFitFASHOS <- list()

for (usePixelInfo in c(FALSE, TRUE)) {
  for (r in 1:repExp) {
    pixInd <- usePixelInfo + 1
    # sample train and test data examples
    nTextures <- length(unique(textureStats$texture))
    sampleTextures <- sample(unique(textureStats$texture))
    trainTextures <- sampleTextures[1:round(nTextures/2)]
    testTextures <- sampleTextures[round(nTextures/2):nTextures]
    trainData <- dplyr::filter(textureStats, texture %in% trainTextures) %>%
      droplevels(.) %>%
      make_task_textures(., nRep)
    testData <- dplyr::filter(textureStats, texture %in% testTextures) %>%
      droplevels(.) %>%
      make_task_textures(., nRep)

    # calculate weights to even out classes
    ratioSame <- sum(trainData$same) / sum(1-trainData$same)
    weights <- rep(1, nrow(trainData))
    weights[trainData$same == 1] <- (1/ratioSame)

    # extract labels
    trainLabel <- trainData$same
    testLabel <- testData$same

    # if selected, do pixel stats model
    if (usePixelInfo) {
      trainStatsPix <- dplyr::select(trainData, all_of(namesPix))
      testStatsPix <- dplyr::select(testData, all_of(namesPix))
      modelFitPix <- cv.glmnet(x = as.matrix(trainStatsPix),
                           y = trainLabel,
                           weights = weights,
                           family="binomial",
                           alpha = alpha,
                           standardize=TRUE,
                           type.measure="class")
      modelPredictionsPix <- predict(modelFitPix, as.matrix(testStatsPix),
                                     type = "class")
      predictionOutcomePix <- mean(as.integer(modelPredictionsPix == testLabel))
    } else {
      modelFitPix <- NA
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
    modelFitFAS <- cv.glmnet(x = as.matrix(trainStatsFAS),
                         y = trainLabel,
                         weights = weights,
                         family="binomial",
                         alpha = alpha,
                         standardize=TRUE,
                         type.measure="class")
    modelPredictionsFAS <- predict(modelFitFAS, as.matrix(testStatsFAS),
                                   type = "class")
    predictionOutcomeFAS <- mean(as.integer(modelPredictionsFAS == testLabel))

    # HOS stats model
    modelFitHOS <- cv.glmnet(x = as.matrix(trainStatsHOS),
                         y = trainLabel,
                         weights = weights,
                         family="binomial",
                         alpha = alpha,
                         standardize=TRUE,
                         type.measure="class")
    modelPredictionsHOS <- predict(modelFitHOS, as.matrix(testStatsHOS), type = "class")
    predictionOutcomeHOS <- mean(as.integer(modelPredictionsHOS == testLabel))

    # all stats model
    modelFitFASHOS <- cv.glmnet(x = as.matrix(trainStatsFASHOS),
                         y = trainLabel,
                         weights = weights,
                         family="binomial",
                         alpha = alpha,
                         standardize=TRUE,
                         type.measure="class")
    modelPredictionsFASHOS <- predict(modelFitFASHOS, as.matrix(testStatsFASHOS), type = "class")
    predictionOutcomeFASHOS <- mean(as.integer(modelPredictionsFASHOS == testLabel))

    results <- c(usePixelInfo, predictionOutcomePix,
                 predictionOutcomeFAS, predictionOutcomeHOS,
                 predictionOutcomeFASHOS)
    resultsDf[nrow(resultsDf)+1,] <- results
  }
}

saveRDS(resultsDf, saveResults)

