library(tidyr)
library(dplyr)
library(glmnet)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_pixNorm.csv"
dataFile <- "../../data/texture_stats/texture_stats_statsNorm.csv"

saveResults <- "../../data/texture_results/texture_dnn.Rds"

usePixelInfo <- TRUE
repExp <- 20
nRep <- 5
layerUnits <- c(30, 10, 3)

# load data
textureStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.)

# permute vector not allowing any value to remain in same place
derangement <- function(x, n){
  outList <- list()
  outList[[1]] <- x
  for (vecN in c(1:n)) {
    incomplete <- TRUE
    while (incomplete) {
      outVec <- base::sample(x)
      allOkTemp <- TRUE
      for (compN in c(1:vecN)) {
        allOkTemp <- allOkTemp * !any(outVec == outList[[compN]])
      }
      incomplete <- incomplete * !allOkTemp
    }
    outList[[vecN+1]] <- outVec
  }
  outList <- outList[-1]
  return(outList)
}


# make a dataframe with substracted statistics of pairs of
# textures that are same or different. For the same pairs,
# it uses all possible combinations of the "quadrant" values.
# For the different pairs, it makes nRep * 2 pairs for each texture
# (that is, nRep pairs in total).
make_task <- function(inputData, nRep){
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
    diffStats <- (reducedDf - scrambledDf) %>%
      dplyr::select(., -texture, - quadrant) %>%
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

# apply the same normalization to trainStats and testStats
normalize_data <- function(trainStats, testStats) {
  trainMat <- as.matrix(trainStats)
  testMat <- as.matrix(testStats)
  columnMeans <- colMeans(trainMat)
  columnSds <- matrixStats::colSds(trainMat)
  if (any(columnSds == 0)) {
    columnSds[which(columnSds == 0)] <- 1
  }
  trainMat <- sweep(trainMat, 2, columnMeans, FUN = "-")
  trainMat <- sweep(trainMat, 2, columnSds, FUN = "/")
  testMat <- sweep(testMat, 2, columnMeans, FUN = "-")
  testMat <- sweep(testMat, 2, columnSds, FUN = "/")
  outStats <- list(trainStats = trainMat, testStats = testMat)
  return(outStats)
}

# do PCA on the stats. Keep 90% variance
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



# run dnn on data
dnn_fit <- function(trainStats, trainLabel, testStats, testLabel, layerUnits) {
  # normalize the input data
  normalizedStats <- normalize_data(trainStats, testStats)
  trainStats <- normalizedStats$trainStats
  testStats <- normalizedStats$testStats
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
parNames <- names(textureStats)
indPix <- c(grep("pix", parNames), grep("^LP", parNames))
namesPix <- parNames[indPix]
indFAS <- c(grep("mm", parNames), grep("acr", parNames))
namesFAS <- parNames[indFAS]
indsHOS <- which(!(parNames %in% c(namesPix, namesFAS)))
indsHOS <- indsHOS[-c(1:2)]
namesHOS <- parNames[indsHOS]

# generate lists and dfs
resultsDf <- data.frame(usePix = logical(), Pix = double(),
                        FA = double(), HOS = double(), FA_HOS = double())
modelFitPix <- list()
modelFitFAS <- list()
modelFitHOS <- list()
modelFitFASHOS <- list()
coefPix <- NULL
coefFAS <- list(NULL, NULL)
coefHOS <- list(NULL, NULL)
coefFASHOS <- list(NULL, NULL)

for (usePixelInfo in c(TRUE, FALSE)) {
  for (r in 1:repExp) {
    pixInd <- usePixelInfo + 1
    nTextures <- length(unique(textureStats$texture))
    sampleTextures <- sample(unique(textureStats$texture))
    trainTextures <- sampleTextures[1:round(nTextures/2)]
    testTextures <- sampleTextures[round(nTextures/2):nTextures]
    trainData <- dplyr::filter(textureStats, texture %in% trainTextures) %>%
      droplevels(.) %>%
      make_task(., nRep)
    testData <- dplyr::filter(textureStats, texture %in% testTextures) %>%
      droplevels(.) %>%
      make_task(., nRep)

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
      predictionOutcomePix <- dnn_fit(trainStatsPix, trainLabel, testStatsPix,
                                      testLabel, layerUnits) 
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
    predictionOutcomeFAS <- dnn_fit(trainStatsFAS, trainLabel, testStatsFAS,
                                    testLabel, layerUnits)
    predictionOutcomeHOS <- dnn_fit(trainStatsHOS, trainLabel, testStatsHOS,
                                    testLabel, layerUnits)
    predictionOutcomeFASHOS <- dnn_fit(trainStatsFASHOS, trainLabel,
                                       testStatsFASHOS, testLabel, layerUnits)

    results <- c(usePixelInfo, predictionOutcomePix, predictionOutcomeFAS,
                 predictionOutcomeHOS, predictionOutcomeFASHOS)
    resultsDf[nrow(resultsDf)+1,] <- results
  }
}

# save results
saveRDS(resultsDf, saveResults)

