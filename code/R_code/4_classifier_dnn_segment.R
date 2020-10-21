library(tidyr)
library(dplyr)
library(glmnet)
library(matrixStats)
source("./analysis_functions.R")
set.seed(2691)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_pixNorm.csv"
dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
saveResults <- "../../data/BSD_results/BSD_dnn2.Rds"

repExp <- 20
layerUnits <- c(30, 10)

# load data
segmentStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.)

# run dnn on data
dnn_fit <- function(trainStats, trainLabel, testStats, testLabel, layerUnits) {
  # normalize the input data
  normalizedStats <- normalize_data(trainStats, testStats)
  trainStats <- normalizedStats$train
  testStats <- normalizedStats$test
  # do PCA on the stats
  statsPCA <- statsPCA(trainStats, testStats, varianceRetained = 0.95)
  trainStats <- statsPCA$trainPCA
  testStats <- statsPCA$testPCA
  # define network
  mod <- kerasR::Sequential()
  mod$add(kerasR::Dense(units = layerUnits[1], input_shape = ncol(trainStats)))
  mod$add(kerasR::ActivityRegularization(l1=0.002))
  mod$add(kerasR::Activation("relu"))
  layerUnits <- layerUnits[-1]
  if (length(layerUnits) == 0) {
    for (u in layerUnits) {
      mod$add(kerasR::Dense(units = u))
      mod$add(kerasR::ActivityRegularization(l1=0.002))
      mod$add(kerasR::Activation("relu"))
    }
  }
  mod$add(kerasR::Dense(units = 1))
  mod$add(kerasR::Activation("sigmoid"))
  kerasR::keras_compile(mod, loss = "binary_crossentropy",
                        metrics = "binary_accuracy", optimizer = kerasR::Adam())
  # calculate weights to even out classes
  ratioSame <- sum(trainLabel) / sum(1-trainLabel)
  weights <- rep(1, nrow(trainData))
  weights[trainLabel == 1] <- (1/ratioSame)
  # fir the model
  kerasR::keras_fit(mod, trainStats, trainLabel,
          batch_size = 32, epochs = 200,
          verbose = 2, validation_data = list(testStats, testLabel))
  # predict classes and return performance
  pred <- kerasR::keras_predict_classes(mod, testStats)
  predictionOutcome <- mean(as.integer(pred == testLabel))
  return(predictionOutcome)
}

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

taskData <- make_task_BSD(segmentStats)
nSegments <- length(unique(segmentStats$ImageName))
imageNames <- unique(segmentStats$ImageName)

for (usePixelInfo in c(TRUE, FALSE)) {
  for (r in 1:repExp) {
    sampleSegments <- sample(imageNames)
    trainSegments <- sampleSegments[1:floor(nSegments*4/5)]
    testSegments <- sampleSegments[ceiling(nSegments*4/5):nSegments]

    trainData <- dplyr::filter(taskData, ImageName %in% trainSegments) %>%
      droplevels(.)
    testData <- dplyr::filter(taskData, ImageName %in% testSegments) %>%
      droplevels(.)

    # extract labels
    trainLabel <- trainData$same
    testLabel <- testData$same

    # if selected, do pixel stats model
    if (usePixelInfo) {
      trainStatsPix <- dplyr::select(trainData, all_of(namesPix))
      testStatsPix <- dplyr::select(testData, all_of(namesPix))
      pixelOutcome <- dnn_fit(trainStatsPix, trainLabel, testStatsPix,
                              testLabel, layerUnits)
    } else {
      modelFitPix <- NA
      predictionOutcomePix <- NA
    }

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

    FASOutcome <- dnn_fit(trainStatsFAS, trainLabel, testStatsFAS,
                          testLabel, layerUnits)
    HOSOutcome <- dnn_fit(trainStatsHOS, trainLabel, testStatsHOS,
                          testLabel, layerUnits)
    FASHOSOutcome <- dnn_fit(trainStatsFASHOS, trainLabel, testStatsFASHOS,
                             testLabel, layerUnits)

    results <- c(usePixelInfo, pixelOutcome, FASOutcome, HOSOutcome, FASHOSOutcome)
    resultsDf[nrow(resultsDf)+1,] <- results
  }
}

# save results
saveRDS(resultsDf, saveResults)

