####################################################
####################################################
#### Functions for data processing and analysis ####
####################################################
####################################################

library(glmnet)
library(dplyr)
library(kerasR)
library(caret)
source("utility_functions.R")

#################################################
###### Functions to operate on the stats matrices
#################################################


# Normalize trainStats and teststats using mean
# and SD from trainstats. If testStats = NA, then only
# normalize and return trainStats
normalize_data <- function(trainStats, testStats=NA) {
  #get means and SD
  trainMeans <- as.numeric(lapply(trainStats, mean))
  trainSD <- as.numeric(lapply(trainStats, sd))
  #scale the stats
  trainStats <- scale(trainStats, center=trainMeans, scale=trainSD)
  if (!is.na(testStats)) {
    testStats <- scale(testStats, center=trainMeans, scale=trainSD)
    normalizedData <- list("train" = trainStats, "test" = testStats)
  } else {
    normalizedData <- trainStats
  }
  return(normalizedData)
}


# do PCA on the stats.
statsPCA <- function(trainStats, testStats, varianceRetained = 0.95) {
  trainPCA <- prcomp(trainStats)
  prcSumm <- summary(trainPCA)
  cummulativeVariance <- prcSumm[["importance"]][3,]
  lastComp <- which(cummulativeVariance >= varianceRetained)[1]
  rotationMat <- prcSumm$rotation[,1:lastComp]
  trainPCAvals <- trainStats %*% rotationMat
  testPCAvals <- testStats %*% rotationMat
  output <- list(trainPCA = trainPCAvals, testPCA = testPCAvals)
  return(output)
}


# do PCA sepparately on the different subgroups of stats
subset_statsPCA <- function(trainStats, testStats, subsetList,
                            varianceRetained=0.95) {
  trainPCall <- data.frame(row.names=c(1:nrow(trainStats)))
  testPCall <- data.frame(row.names=c(1:nrow(testStats)))
  if (!is.list(subsetList)) {
    subsetList <- list(all=subsetList)
  }
  subsetNames <- names(subsetList) 
  for (t in c(1:length(subsetList))) {
    t_name <- subsetNames[t]
    subsetStats <- subsetList[[t]]
    # check that the subset list stats are in the train stats
    if (any(subsetStats %in% names(trainStats))) {
      subsetTrain <- dplyr::select(trainStats, all_of(subsetStats))
      subsetTest <- dplyr::select(testStats, all_of(subsetStats))
      trainPCA <- prcomp(subsetTrain)
      prcSumm <- summary(trainPCA)
      cummulativeVariance <- prcSumm[["importance"]][3,]
      lastComp <- which(cummulativeVariance >= varianceRetained)[1]
      rotationMat <- prcSumm$rotation[,1:lastComp]
      trainPCAvals <- as.data.frame(as.matrix(subsetTrain) %*% rotationMat)
      testPCAvals <- as.data.frame(as.matrix(subsetTest) %*% rotationMat)
      # put the PCA vals into a dataframe
      variableNames <- paste(t_name, "_PC", c(1:lastComp), sep="")
      names(trainPCAvals) <- variableNames
      names(testPCAvals) <- variableNames
      trainPCall <- cbind(trainPCall, trainPCAvals)
      testPCall <- cbind(testPCall, testPCAvals)
    }
  }
  output <- list(trainPCA = trainPCall, testPCA = testPCall)
  return(output)
}


#################################################
###### Functions to preprocess data for analysis
#################################################

# Make the dataframe with the task for the BSD dataset.
# Generates 2 different kinds of pairs:
# Pairs from different segments (the center segments paired
# with each neighbor), and pairs from same segment (a given
# segment is split in half, be it center segment or a neighbor)
make_task_BSD <- function(inputData){
  # make the rows for the "different" condition
  diffDf <- inputData %>%
   dplyr::filter(., Type %in% c("center", "neighbor")) %>%
   group_by(., ImageName) %>%
   dplyr::mutate_each(funs(abs(. - .[Type == "center"])),
                      -Type, -Segment) %>%
   dplyr::filter(., Type != "center") %>%
   ungroup(.) %>%
   cbind(same = 0, .)
  # make the rows for the "same" condition
  sameDfCenter <- inputData %>%
   dplyr::filter(., Type == "half_center") %>%
   group_by(., ImageName) %>%
   dplyr::mutate_each(funs(abs(. - .[Segment == 1])),
                      -Type, -Segment) %>%
   dplyr::filter(., Segment != 1) %>%
   ungroup(.) %>%
   cbind(same = 1, .)
  sameDfNeighbor <- inputData %>%
   dplyr::filter(., Type == "half_neighbor") %>%
   group_by(., ImageName, Segment) %>%
   dplyr::mutate(., half = c(1:2)) %>%
   dplyr::mutate_each(funs(abs(. - .[half == 1])),
                      -Type, -half) %>%
   dplyr::filter(., half != 1) %>%
   dplyr::select(., -half) %>%
   ungroup(.) %>%
   cbind(same = 1, .)
  allData <- rbind(sameDfCenter, sameDfNeighbor, diffDf)
  return(allData)
}


# Make the dataframe with the task for the textures.
# For the same pairs,
# it uses all possible combinations of the "quadrant" values.
# For the different pairs, it makes nRep * 2 pairs for each texture
# (that is, nRep pairs in total).
make_task_textures <- function(inputData, nRep){
  # make the rows for the "same" condition
  inputQuadrants <- unique(inputData$quadrant)
  quadrantsCombinations <- combn(inputQuadrants, 2)
  splitData <- dplyr::select(inputData, -texture, -quadrant) %>%
    split(., inputData$quadrant)
  sameDf <- data.frame()
  for (col in c(1:ncol(quadrantsCombinations))) {
    quadrants <- quadrantsCombinations[,col]
    sameStats <- as_tibble(splitData[[as.character(quadrants[1])]] -
                           splitData[[as.character(quadrants[2])]]) %>%
      abs(.)
    sameCond <- data.frame(texture1 = unique(inputData$texture),
                           texture2 = unique(inputData$texture), same = 1)
    sameDf <- rbind(sameDf, cbind(sameCond, sameStats)) %>%
      as_tibble(.)
  }
  # combine different textures for the "different" condition
  # sample different combinations that are not repeated
  diffSorts <- derangement(c(1:length(unique(inputData$texture))), nRep)
  diffDf <- data.frame()
  for (n in c(1:nRep)) {
    sampleQuadrant <- sample(unique(inputData$quadrant), 1)
    # substract the statistics
    reducedDf <- dplyr::filter(inputData, quadrant == sampleQuadrant)
    scrambledDf <- reducedDf[diffSorts[[n]],]  
    # put statistics and condition info into a dataframe
    reducedStats <- dplyr::select(reducedDf, -texture, -quadrant)
    scrambledStats <- dplyr::select(scrambledDf, -texture, -quadrant)
    diffStats <- (reducedStats - scrambledStats) %>%
      abs(.)
    diffCond <- data.frame(texture1 = unique(reducedDf$texture),
                           texture2 = unique(scrambledDf$texture), same = 0)
    diffDf <- rbind(diffDf, cbind(diffCond, diffStats))
  }
  # put same and different dataframes together   
  allData <- rbind(sameDf, diffDf)
  allData$same <- as.integer(allData$same)
  return(allData)
}


###############################
###### Analysis functionns
###############################

# make a dataframe with dot products of pairs of
# textures that belong to the same image in inputData.
# Recofnizes Type and compares center to neighbors and the
# two halves of the center.
compute_angles <- function(inputData){
  imageNames <- unique(inputData$ImageName)
  angleDf <- data.frame(ImageName = character(), Type = character(),
                        Segment = numeric(), angle = numeric())
  for (im in c(1:length(imageNames))) {
    imageData <- dplyr::filter(inputData, ImageName == imageNames[im])
    # angle between neighbors
    centerVec <- dplyr::filter(imageData, Type == "center") %>%
      dplyr::select(., -ImageName, -Type, -Segment) %>%
      as.matrix(.)
    neighborVec <- dplyr::filter(imageData, Type == "neighbor") %>%
      dplyr::select(., -ImageName, -Type, -Segment) %>%
      as.matrix(.)
    neighborAngle <- vec_angle(centerVec, neighborVec)
    # angle between center halves
    splitCenterVec <- dplyr::filter(imageData, Type == "half_center") %>%
      dplyr::select(., -ImageName, -Type, -Segment) %>%
      as.matrix(.)
    centerAngle <- vec_angle(splitCenterVec[1,], splitCenterVec[2,])
    # angle between center halves
    halfNeighborAngle <- NULL
    splitNeighborDf <- dplyr::filter(imageData, Type == "half_neighbor") 
    for (seg in unique(splitNeighborDf$Segment)) {
      splitNeighborVec <- dplyr::filter(splitNeighborDf, Segment == seg) %>%
        dplyr::select(., -ImageName, -Type, -Segment) %>%
        as.matrix(.)
      halfNeighborAngle <- c(halfNeighborAngle, vec_angle(splitNeighborVec[1,],
                                                              splitNeighborVec[2,]))
    }
    # make data frame
    nNeighbors <- length(neighborAngle)
    nHalfNeighbors <- length(halfNeighborAngle)
    nameVec <- rep(imageNames[im], nNeighbors + nHalfNeighbors + 1)
    typeVec <- c(rep("diff", nNeighbors), rep("same", 1 + nHalfNeighbors))
    segmentVec <- c(1:(nNeighbors + nHalfNeighbors), 0)
    angles <- c(neighborAngle, halfNeighborAngle, centerAngle)
    tempDf <- data.frame(ImageName = nameVec, Type = typeVec,
                         Segment = segmentVec, angle = angles)
    angleDf <- rbind(angleDf, tempDf)
  }
  return(angleDf)
}

# Generate a list with the prepared data to train and test a model
prepare_data_fit_test <- function(trainData, testData, statsToUse=NA,
                             balanceWeights=TRUE, subsetsPCA=NA,
                             labelColumn="same") {
  # calculate weights to even out classes
  weights <- rep(1, nrow(trainData))
  if (balanceWeights) {
    ratioSame <- sum(trainData[[labelColumn]]) / sum(1-trainData[[labelColumn]])
    weights[trainData[[labelColumn]] == 1] <- (1/ratioSame)
  }
  # extract labels
  trainLabel <- trainData[[labelColumn]]
  testLabel <- testData[[labelColumn]]
  # extract statistics
  if (is.na(statsToUse[1])) {
    allStatsNames <- get_statistics_names(names(trainData))
    statsToUse <- c(allStatsNames$pixel, allStatsNames$FAS, allStatsNames$HOS)  
  } else if (is.list(statsToUse)) {
      statsToUse <- unlist_names(statsToUse)
  }
  trainStats <- dplyr::select(trainData, all_of(statsToUse))
  testStats <- dplyr::select(testData, all_of(statsToUse))
  # Normalize the statistics
  normalizedStats <- normalize_data(trainStats=trainStats, testStats=testStats)
  trainStats <- as.data.frame(normalizedStats$train)
  testStats <- as.data.frame(normalizedStats$test)
  # If required, do PCA
  if (!is.na(subsetsPCA[1])) {
    pcaStats <- subset_statsPCA(trainStats, testStats, subsetsPCA)
    trainStats <- pcaStats$trainPCA
    testStats <- pcaStats$testPCA
  }
  # fit the model
  preparedData <- list(trainStats=trainStats,
                       testStats=testStats,
                       trainLabel=trainLabel,
                       testLabel=testLabel,
                       weights=weights)
  return(preparedData)
}


# Train a ridge regression model on given train and test data.
train_test_ridge <- function(trainData, testData, statsToUse=NA,
                             balanceWeights = TRUE, subsetsPCA=NA,
                             labelColumn="same") {
  preparedData <- prepare_data_fit_test(trainData=trainData,
                                        testData=testData,
                                        statsToUse=statsToUse,
                                        balanceWeights=balanceWeights,
                                        subsetsPCA=subsetsPCA,
                                        labelColumn=labelColumn)
  # fit the model
  modelFit <- cv.glmnet(x = as.matrix(preparedData$trainStats),
                       y = preparedData$trainLabel,
                       weights = preparedData$weights,
                       family="binomial",
                       alpha = 0, # alpha = 0 is ridge regression
                       standardize=FALSE,
                       type.measure="class")
  modelPredictions <- predict(modelFit, as.matrix(preparedData$testStats),
                                 type = "class")
  modelResponses <- predict(modelFit, as.matrix(preparedData$testStats),
                                 type = "link")
  correctPredictions <- as.integer(modelPredictions == preparedData$testLabel)
  predictionOutcome <- mean(correctPredictions)
  referenceLabels <- factor(preparedData$testLabel, levels=c("0","1"))
  confusionMatrix <- caret::confusionMatrix(factor(modelPredictions),
                                                     referenceLabels)
  modelOutput <- list(predictions = modelPredictions,
                      accuracy = predictionOutcome,
                      correctPredictions = correctPredictions,
                      confusionMatrix = confusionMatrix,
                      modelResponses = modelResponses)
  return(modelOutput)
}


# make dnn model
# tutorial: https://cran.r-project.org/web/packages/kerasR/vignettes/introduction.html
make_dnn_model <- function(layerUnits, inputShape, regularizationWeight) {
  mod <- kerasR::Sequential()
  # only the first layer requires input_shape
  mod$add(kerasR::Dense(units=layerUnits[1], input_shape=inputShape))
  mod$add(kerasR::ActivityRegularization(l1=regularizationWeight))
  mod$add(kerasR::Activation("relu"))
  layerUnits <- layerUnits[-1]
  if (length(layerUnits)!=0) {
    for (u in layerUnits) {
      mod$add(kerasR::Dense(units=u))
      mod$add(kerasR::ActivityRegularization(l1=regularizationWeight))
      mod$add(kerasR::Activation("relu"))
    }
  }
  mod$add(kerasR::Dense(units=1))
  mod$add(kerasR::Activation("sigmoid"))
  kerasR::keras_compile(mod, loss="binary_crossentropy",
                        metrics="binary_accuracy", optimizer=kerasR::Adam())
  return(mod)
}


# Train a ridge regression model on given train and test
# data.
train_test_dnn <- function(trainData, testData, statsToUse=NA,
                           balanceWeights=TRUE, subsetsPCA=NA,
                           layerUnits=c(30,10), regularizationWeight=0.002,
                           epochs=200) {
  preparedData <- prepare_data_fit_test(trainData, testData,
                                        statsToUse, balanceWeights,
                                        subsetsPCA)
  modelDNN <- make_dnn_model(layerUnits, ncol(preparedData$trainStats),
                             regularizationWeight)
  # move sample weights to python
  np <- reticulate::import("numpy", convert=FALSE)
  pyWeights <- np$array(reticulate::r_to_py(preparedData$weights))
  # fit the model
  kerasR::keras_fit(modelDNN, as.matrix(preparedData$trainStats),
                    preparedData$trainLabel, batch_size=32,
                    epochs = epochs, verbose = 2,
                    validation_data = list(as.matrix(preparedData$testStats),
                                           preparedData$testLabel),
                    sample_weight=pyWeights)
                    #class_weight=list(0.1,1))
  # predict classes and return performance
  modelPredictions <- kerasR::keras_predict_classes(modelDNN,
                                        as.matrix(preparedData$testStats))
  correctPredictions <- as.integer(modelPredictions == preparedData$testLabel)
  predictionOutcome <- mean(correctPredictions)
  confusionMatrix <- caret::confusionMatrix(factor(modelPredictions),
                                                     factor(preparedData$testLabel))
  modelOutput <- list(predictions=modelPredictions,
                      accuracy=predictionOutcome,
                      correctPredictions=correctPredictions,
                      confusionMatrix=confusionMatrix,
                      accuracyHistory=modelDNN$history$history$val_binary_accuracy)
  return(modelOutput)
}

