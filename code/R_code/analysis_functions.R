#### Functions for data processing and analysis ####

#### Multi-use

# from a vector with names of columns of a dataframe containing
# Portilla-Simoncelli statistics, separate the names into the
# kinds of statistics and experiment design names
get_statistics_names <- function(dfNamesVec) {
  # get the indices of the different types of stats to use
  indPix <- c(grep("pix", dfNamesVec), grep("^LP", dfNamesVec))
  namesPix <- dfNamesVec[indPix]
  indFAS <- c(grep("mm", dfNamesVec), grep("acr", dfNamesVec))
  namesFAS <- dfNamesVec[indFAS]
  indHOS <- c(grep("acm", dfNamesVec), grep("cmc", dfNamesVec),
               grep("pmc", dfNamesVec), grep("prc", dfNamesVec))
  namesHOS <- dfNamesVec[indHOS]
  indDesign <- which(!(dfNamesVec %in% c(namesPix, namesFAS, namesHOS)))
  namesDesign <- dfNamesVec[indDesign]
  subsetNames <- list(pixel = namesPix, FAS = namesFAS, HOS = namesHOS,
                      design = namesDesign)
  return(subsetNames)
}


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

# normalize data so that trainStats have SD = 1
# and mean 0, and testStats are normalized with the
# same parameters. If testStats = NA, then only normalize
# the trainStats
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

# make the task for the BSD dataset generating the pairs
# with the differences in statistics. Generates 3 different
# kinds of pairs: Within the pairs from the same segment
# there is the split center segment and the split neighboring
# segments, and the the pairs of different segments (center and
# neighbors)
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


# make a dataframe with substracted statistics of pairs of
# textures that are same or different. For the same pairs,
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

# do PCA sepparately on the different subgroups of stats
subset_statsPCA <- function(trainStats, testStats, varianceRetained = 0.95) {
  # get the subsets to which each stat belongs 
  statsNames <- names(trainStats)
  statsType <- subset_statistics(statsNames = statsNames)
  types <- unique(statsType)
  # get design cols
  namesTemp <- get_statistics_names(statsNames)
  designNames <- namesTemp$design
  # get the PCs for each subset separtely
  trainPCall <- data.frame(row.names = c(1:nrow(trainStats)))
  testPCall <- data.frame(row.names = c(1:nrow(testStats)))
  for (t in types) {
    subsetStats <- statsNames[which(statsType == t)]
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
    variableNames <- paste(t, "_PC", c(1:lastComp), sep="")
    names(trainPCAvals) <- variableNames
    names(testPCAvals) <- variableNames
    trainPCall <- cbind(trainPCall, trainPCAvals)
    testPCall <- cbind(testPCall, testPCAvals)
  }
  output <- list(trainPCA = trainPCall, testPCA = testPCall)
  return(output)
}

# takes a vector of class predictions and a reference vector, and
# returns a vector indicating for each prediction whether it is
# a true positive (TP), false negative (FN), true negative (TN),
# or false positive (FP)
confusion_indices <- function(predictions, reference) {
  trueInds <- predictions == reference
  falseInds <- predictions != reference
  positiveInd <- as.integer(reference) == 1
  negativeInd <- as.integer(reference) != 1
  confusionVec <- NULL
  confusionVec[which(trueInds & positiveInd)] <- "TP"
  confusionVec[which(falseInds & positiveInd)] <- "FN"
  confusionVec[which(trueInds & negativeInd)] <- "TN"
  confusionVec[which(falseInds & negativeInd)] <- "FP"
  return(confusionVec)
}

# return a vector indicating the subset to which each statistic
# in a statistics name vector belongs
subset_statistics <- function(statsNames) {
  paramTypeNames <- c("mm", "acr", "acm", "cmc", "pmc", "prc")
  tempVec <- rep(NA, length(statsNames))
  for (pN in paramTypeNames) {
    paramInds  <- grep(paste("^", pN, sep=""), statsNames)
    tempVec[paramInds] <- pN
  }
  return(tempVec)
}


####################
###### used for ROC
####################

# calculate angle between vectors. The second input
# can be a matrix, and the dot product is computed
# with each row
vec_angle <- function(inputVec, inputMat) {
  if (is.null(dim(inputMat))) {
    dotProd <-  inputVec %*% inputMat
    matLength <- sqrt(sum(inputMat^2))
  } else {
    dotProd <-  inputVec %*% t(inputMat)
    matLength <- sqrt(rowSums(inputMat^2))
  }
  vecLength <- sqrt(sum(inputVec^2))
  angleCos <- dotProd / (vecLength * matLength)
  angles <- acos(pmin(pmax(angleCos, -1.0), 1.0))
  return(as.numeric(angles))
}
  

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

