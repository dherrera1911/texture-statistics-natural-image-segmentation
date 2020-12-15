######################################################
######################################################
#### Functions providing simple functionanilities ####
######################################################
######################################################


# separate the columns names of a dataframe between
# the design columns and the different kinds of
# statistics. If subsetHOS = TRUE, in the HOS field
# of the output returns a list with all subsets of
# HOS
get_statistics_names <- function(dfNamesVec, subsetHOS = FALSE) {
  # get the indices of the different types of stats to use
  indPix <- c(grep("pix", dfNamesVec), grep("^LP", dfNamesVec))
  namesPix <- dfNamesVec[indPix]
  indFAS <- c(grep("mm", dfNamesVec), grep("acr", dfNamesVec))
  namesFAS <- dfNamesVec[indFAS]
  if (!subsetHOS) {
    indHOS <- c(grep("acm", dfNamesVec), grep("cmc", dfNamesVec),
                 grep("pmc", dfNamesVec), grep("prc", dfNamesVec))
    namesHOS <- dfNamesVec[indHOS]
    namesHOSVec <- namesHOS
  } else {
    indACM <- grep("acm", dfNamesVec)
    indCMC <- grep("cmc", dfNamesVec)
    indPMC <- grep("pmc", dfNamesVec)
    indPRC <- grep("prc", dfNamesVec)
    namesHOS <- list(acm = dfNamesVec[indACM],
                     cmc = dfNamesVec[indCMC],
                     pmc = dfNamesVec[indPMC],
                     prc = dfNamesVec[indPRC])
    namesHOSVec <- unlist_names(namesHOS)
  }
  indDesign <- which(!(dfNamesVec %in% c(namesPix, namesFAS, namesHOSVec)))
  namesDesign <- dfNamesVec[indDesign]
  subsetNames <- list(pixel = namesPix, FAS = namesFAS, HOS = namesHOS,
                      design = namesDesign)
  return(subsetNames)
}

# remove statistics with 0 variance
remove_constant_stats <- function(rawData) {
  varNames <- get_statistics_names(names(rawData))
  statisticsNames <- varNames[which(names(varNames)!="design")]
  statisticsNames <- unlist_names(statisticsNames)
  statsDf <- dplyr::select(rawData, statisticsNames) 
  statsSD <- as.numeric(lapply(statsDf, sd))
  noVarInds <- which(statsSD == 0)
  if (length(noVarInds)>0) {
    statsDf <- statsDf[, -noVarInds]
  }
  designDf <- dplyr::select(rawData, varNames$design)
  outputDf <- cbind(designDf, statsDf)
  return(outputDf)
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


unlist_names <- function(namesList) {
  namesVec <- NULL
  for (l in c(1:length(namesList))) {
    namesVec <- c(namesVec, namesList[[l]])
  }
  return(namesVec)
}

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


