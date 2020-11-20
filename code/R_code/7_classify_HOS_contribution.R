library(tidyr)
library(dplyr)
library(glmnet)
library(caret)
library(vcd)
source("./analysis_functions.R")
set.seed(2691)

alpha <- 0
dataSplits <- 10
doPCA <- TRUE

# data files
dataFileHOS <- "../../data/BSD_stats/pairs_HOS_better.Rds"
dataFileFASHOS <- "../../data/BSD_stats/pairs_FASHOS_better.Rds"
dataFileSegments <- "../../data/BSD_stats/BSD_stats_Corr.csv"
# saving files
savingFileHOS <- "../../data/BSD_stats/analysis_better_HOS_PCA.RDS"
savingFileFASHOS <- "../../data/BSD_stats/analysis_better_FASHOS_PCA.RDS"
savingFileDiff <- "../../data/BSD_stats/analysis_better_FASHOS_Diff_PCA.RDS"

betterHOSDf <- readRDS(dataFileHOS)
betterFASHOSDf <- readRDS(dataFileFASHOS)
segmentStats <- read.csv(dataFileSegments, sep=",") %>%
  as_tibble(.)

# get the indices of the different types of stats to use
parNames <- names(segmentStats)
statisticsNames <- get_statistics_names(parNames)
namesPix <- statisticsNames$pixel
namesFAS <- statisticsNames$FAS
namesHOS <- statisticsNames$HOS
designNames <- statisticsNames$design

################
# Do linear regression to see whether the pairs with better HOS are clearly
# separated
################

testSetSplits <- 10
allData <- make_task_BSD(segmentStats) %>%
  as_tibble(.)

# Make a dataframe with all pairs of segments indicating if
# HOS performed better. Separate between test and training set
# pairs of images where HOS was not better
nonBetterHOS <- anti_join(allData, betterHOSDf, by = c("same", "ImageName", "Type", "Segment"))
# add labels
nonBetterHOS$betterHOS <- 0
betterHOSDf$betterHOS <- 1

# make splits of images to test for better HOS
nBetterHOSSegments <- nrow(betterHOSDf)
splitSize <- ceiling(nBetterHOSSegments/testSetSplits)
segmentListHOS <- split(sample(nBetterHOSSegments),
                          ceiling(seq_along(c(1:nBetterHOSSegments))/splitSize))
nNonHOSSegments <- nrow(nonBetterHOS)
splitSize <- ceiling(nNonHOSSegments/testSetSplits)
segmentListNonHOS <- split(sample(nNonHOSSegments),
                          ceiling(seq_along(c(1:nNonHOSSegments))/splitSize))

coefBetterHOS <- list()
confusionMatrixHOS <- list()
confusionHOSDf <- NULL

for (sp in c(1:testSetSplits)) {
  # separate HOS Df into train and test
  testHOS <- betterHOSDf[segmentListHOS[[sp]],]
  trainHOS <- betterHOSDf[-segmentListHOS[[sp]],]
  # separate FAS Df into train and test
  testFAS <- nonBetterHOS[segmentListNonHOS[[sp]],]
  trainFAS <- nonBetterHOS[-segmentListNonHOS[[sp]],]
  # bind into complete train and test Dfs
  trainDf <- rbind(trainHOS, trainFAS)
  testDf <- rbind(testHOS, testFAS)
  # separate stats and labels
  trainStats <- dplyr::select(trainDf, all_of(c(namesFAS, namesHOS)))
  testStats <- dplyr::select(testDf, all_of(c(namesFAS, namesHOS)))
  statsNamesTemp <- names(trainStats)
  normalizedStats <- normalize_data(trainStats = trainStats,
                                          testStats = testStats)
  trainStats <- as.data.frame(normalizedStats$train)
  testStats <- as.data.frame(normalizedStats$test)
  names(trainStats) <- statsNamesTemp
  names(testStats) <- statsNamesTemp
  # if do PCA, convert the stats into PCs for each subset of
  # statistics
  if (doPCA) {
    statsPCA <- subset_statsPCA(trainStats, testStats)
    trainStats <- statsPCA$trainPCA
    testStats <- statsPCA$testPCA
  }

  # get the labels
  trainLabels <- trainDf$betterHOS
  testLabels <- testDf$betterHOS

  # get vector of weights for training
  ratioHOS <- sum(trainLabels)/sum(1-trainLabels)
  weights <- rep(1, length(trainLabels))
  weights[trainLabels==1] <- (1/ratioHOS)

  # train a linear discriminator model
  modelBetterHOS <- cv.glmnet(x = as.matrix(trainStats),
                       y = trainLabels,
                       weights = weights,
                       family="binomial",
                       alpha = alpha,
                       standardize=FALSE,
                       type.measure="class")
  modelPredictionsBetterHOS <- predict(modelBetterHOS,
                                    as.matrix(testStats), type = "class")
  predictionOutcomeBetterHOS <- mean(as.integer(modelPredictionsBetterHOS ==
                                                testLabels))
  coefBetterHOS[[sp]] <- coef(modelBetterHOS)

  # Do a sensitivity specificity confusion matrix analysis
  predictionsFactor <- factor(modelPredictionsBetterHOS, levels = c("1", "0"))
  referenceFactor <- factor(testLabels, levels = c("1", "0"))
  confusionMatrixHOS[[sp]] <- caret::confusionMatrix(predictionsFactor,
                                                     referenceFactor)
  # store in a Df the status of each tested pair
  confusionVec <- confusion_indices(predictionsFactor, referenceFactor)
  confusionHOSDf <- dplyr::select(testDf, same, ImageName, Type, Segment) %>%
    dplyr::mutate(., classificationOutcomeHOS = confusionVec) %>%
    rbind(confusionHOSDf, .)
}


# get the global confusion matrix pooling the models
globalMatrixHOS <- confusionMatrixHOS[[1]]$table
for (sp in c(2:length(confusionMatrixHOS))) {
  globalMatrixHOS <- globalMatrixHOS + confusionMatrixHOS[[sp]]$table
}

kappaHOS <- vcd::Kappa(globalMatrixHOS)

statisticsHOSBetter <- list(kappaHOS = kappaHOS, globalMatrixHOS = globalMatrixHOS,
                            confusionHOSDf = confusionHOSDf, 
                            coefBetterHOS = coefBetterHOS)

saveRDS(statisticsHOSBetter, savingFileHOS)


################
# Do linear regression to see whether the pairs with better FASHOS are clearly
# separated
################

# pairs of images where HOS was not better
nonBetterFASHOS <- anti_join(allData, betterFASHOSDf, by = c("same", "ImageName", "Type", "Segment"))
# add labels
nonBetterFASHOS$betterHOS <- 0
betterFASHOSDf$betterHOS <- 1

# make splits of images to test for better HOS
nBetterFASHOSSegments <- nrow(betterFASHOSDf)
splitSize <- ceiling(nBetterFASHOSSegments/testSetSplits)
segmentListFASHOS <- split(sample(nBetterFASHOSSegments),
                          ceiling(seq_along(c(1:nBetterFASHOSSegments))/splitSize))
nNonFASHOSSegments <- nrow(nonBetterFASHOS)
splitSize <- ceiling(nNonFASHOSSegments/testSetSplits)
segmentListNonFASHOS <- split(sample(nNonFASHOSSegments),
                          ceiling(seq_along(c(1:nNonFASHOSSegments))/splitSize))

coefBetterFASHOS <- list()
confusionMatrixFASHOS <- list()
confusionFASHOSDf <- NULL

for (sp in c(1:testSetSplits)) {
  # separate HOS Df into train and test
  testFASHOS <- betterFASHOSDf[segmentListFASHOS[[sp]],]
  trainFASHOS <- betterFASHOSDf[-segmentListFASHOS[[sp]],]
  # separate FAS Df into train and test
  testFAS <- nonBetterFASHOS[segmentListNonFASHOS[[sp]],]
  trainFAS <- nonBetterFASHOS[-segmentListNonFASHOS[[sp]],]
  # bind into complete train and test Dfs
  trainDf <- rbind(trainFASHOS, trainFAS)
  testDf <- rbind(testFASHOS, testFAS)
  # separate stats and labels
  trainStats <- dplyr::select(trainDf, all_of(c(namesFAS, namesHOS)))
  testStats <- dplyr::select(testDf, all_of(c(namesFAS, namesHOS)))
  statsNamesTemp <- names(trainStats)
  normalizedStats <- normalize_data(trainStats = trainStats,
                                          testStats = testStats)
  trainStats <- as.data.frame(normalizedStats$train)
  testStats <- as.data.frame(normalizedStats$test)
  names(trainStats) <- statsNamesTemp
  names(testStats) <- statsNamesTemp

  # if do PCA, convert the stats into PCs for each subset of
  # statistics
  if (doPCA) {
    statsPCA <- subset_statsPCA(trainStats, testStats)
    trainStats <- statsPCA$trainPCA
    testStats <- statsPCA$testPCA
  }

  # get the labels
  trainLabels <- trainDf$betterHOS
  testLabels <- testDf$betterHOS

  # get vector of weights for training
  ratioFASHOS <- sum(trainLabels)/sum(1-trainLabels)
  weights <- rep(1, length(trainLabels))
  weights[trainLabels==1] <- (1/ratioFASHOS)

  # train a linear discriminator model
  modelBetterFASHOS <- cv.glmnet(x = as.matrix(trainStats),
                       y = trainLabels,
                       weights = weights,
                       family="binomial",
                       alpha = alpha,
                       standardize=FALSE,
                       type.measure="class")
  modelPredictionsBetterFASHOS <- predict(modelBetterFASHOS,
                                    as.matrix(testStats), type = "class")
  predictionOutcomeBetterFASHOS <- mean(as.integer(modelPredictionsBetterFASHOS ==
                                                   testLabels))
  coefBetterFASHOS[[sp]] <- coef(modelBetterFASHOS)

  # Do a sensitivity specificity confusion matrix analysis
  predictionsFactor <- factor(modelPredictionsBetterFASHOS, levels = c("1", "0"))
  referenceFactor <- factor(testLabels, levels = c("1", "0"))
  confusionMatrixFASHOS[[sp]] <- caret::confusionMatrix(predictionsFactor,
                                                        referenceFactor)
  # store in a Df the status of each tested pair
  confusionVec <- confusion_indices(predictionsFactor, referenceFactor)
  confusionFASHOSDf <- dplyr::select(testDf, same, ImageName, Type, Segment) %>%
    dplyr::mutate(., classificationOutcomeFASHOS = confusionVec) %>%
    rbind(confusionFASHOSDf, .)
}

# get the global confusion matrix pooling the models
globalMatrixFASHOS <- confusionMatrixFASHOS[[1]]$table
for (sp in c(2:length(confusionMatrixFASHOS))) {
  globalMatrixFASHOS <- globalMatrixFASHOS + confusionMatrixFASHOS[[sp]]$table
}

kappaFASHOS <- vcd::Kappa(globalMatrixFASHOS)

statisticsFASHOSBetter <- list(kappaFASHOS = kappaFASHOS,
                               globalMatrixFASHOS = globalMatrixFASHOS,
                               confusionFASHOSDf = confusionFASHOSDf,
                               coefBetterFASHOS = coefBetterFASHOS)

saveRDS(statisticsFASHOSBetter, savingFileFASHOS)



##################
## Analyze performance for only pairs that are not same segment
#################

diffPairsDf <- dplyr::filter(allData, same == 0)
betterDiffDf <- dplyr::filter(betterFASHOSDf, same == 0)

# pairs of images where HOS was not better
nonBetterDiff <- anti_join(diffPairsDf, betterDiffDf, by = c("same", "ImageName", "Type", "Segment"))
nonBetterDiff$betterHOS <- 0
betterDiffDf$betterHOS <- 1

# make splits of images to test for better HOS
nBetterDiffSegments <- nrow(betterDiffDf)
splitSize <- ceiling(nBetterDiffSegments/testSetSplits)
segmentListDiff <- split(sample(nBetterDiffSegments),
                          ceiling(seq_along(c(1:nBetterDiffSegments))/splitSize))
nNonDiffSegments <- nrow(nonBetterDiff)
splitSize <- ceiling(nNonDiffSegments/testSetSplits)
segmentListNonDiff <- split(sample(nNonDiffSegments),
                          ceiling(seq_along(c(1:nNonDiffSegments))/splitSize))

coefBetterDiff <- list()
confusionMatrixDiff <- list()
confusionDiffDf <- NULL

for (sp in c(1:testSetSplits)) {
  # separate HOS Df into train and test
  testDiff <- betterDiffDf[segmentListDiff[[sp]],]
  trainDiff <- betterDiffDf[-segmentListDiff[[sp]],]
  # separate FAS Df into train and test
  testFAS <- nonBetterDiff[segmentListNonDiff[[sp]],]
  trainFAS <- nonBetterDiff[-segmentListNonDiff[[sp]],]
  # bind into complete train and test Dfs
  trainDf <- rbind(trainDiff, trainFAS)
  testDf <- rbind(testDiff, testFAS)
  # separate stats and labels
  trainStats <- dplyr::select(trainDf, all_of(c(namesFAS, namesHOS)))
  testStats <- dplyr::select(testDf, all_of(c(namesFAS, namesHOS)))
  statsNamesTemp <- names(trainStats)
  normalizedStats <- normalize_data(trainStats = trainStats,
                                          testStats = testStats)
  trainStats <- as.data.frame(normalizedStats$train)
  testStats <- as.data.frame(normalizedStats$test)
  names(trainStats) <- statsNamesTemp
  names(testStats) <- statsNamesTemp
  # if do PCA, convert the stats into PCs for each subset of
  # statistics
  if (doPCA) {
    statsPCA <- subset_statsPCA(trainStats, testStats)
    trainStats <- statsPCA$trainPCA
    testStats <- statsPCA$testPCA
  }

  # get the labels
  trainLabels <- trainDf$betterHOS
  testLabels <- testDf$betterHOS

  # get vector of weights for training
  ratioDiff <- sum(trainLabels)/sum(1-trainLabels)
  weights <- rep(1, length(trainLabels))
  weights[trainLabels==1] <- (1/ratioDiff)

  # train a linear discriminator model
  modelBetterDiff <- cv.glmnet(x = as.matrix(trainStats),
                       y = trainLabels,
                       weights = weights,
                       family="binomial",
                       alpha = alpha,
                       standardize=FALSE,
                       type.measure="class")
  modelPredictionsBetterDiff <- predict(modelBetterDiff,
                                    as.matrix(testStats), type = "class")
  predictionOutcomeBetterDiff <- mean(as.integer(modelPredictionsBetterDiff == testLabels))
  coefBetterDiff[[sp]] <- coef(modelBetterDiff)

  # Do a sensitivity specificity confusion matrix analysis
  predictionsFactor <- factor(modelPredictionsBetterDiff, levels = c("1", "0"))
  referenceFactor <- factor(testLabels, levels = c("1", "0"))
  confusionMatrixDiff[[sp]] <- caret::confusionMatrix(predictionsFactor,
                                                      referenceFactor)
  # store in a Df the status of each tested pair
  confusionVec <- confusion_indices(predictionsFactor, referenceFactor)
  confusionDiffDf <- dplyr::select(testDf, same, ImageName, Type, Segment) %>%
    dplyr::mutate(., classificationOutcomeDiff = confusionVec) %>%
    rbind(confusionDiffDf, .)
}

# get the global confusion matrix pooling the models
globalMatrixDiff <- confusionMatrixDiff[[1]]$table
for (sp in c(2:length(confusionMatrixDiff))) {
  globalMatrixDiff <- globalMatrixDiff + confusionMatrixDiff[[sp]]$table
}

kappaDiff <- vcd::Kappa(globalMatrixDiff)

statisticsDiffBetter <- list(kappaDiff = kappaDiff,
                             globalMatrixDiff = globalMatrixDiff,
                             confusionDiffDf = confusionDiffDf,
                             coefBetterDiff = coefBetterDiff)

saveRDS(statisticsDiffBetter, savingFileDiff)

