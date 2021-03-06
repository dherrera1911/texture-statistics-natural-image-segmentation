library(tidyr)
library(dplyr)
library(glmnet)
library(caret)
library(vcd)
source("./analysis_functions.R")
set.seed(2691)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_statsNorm.csv"
#saveAgreementFile <- "../../data/BSD_results/6_classification_agreement_texture.Rds"
dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
saveAgreementFile <- "../../data/BSD_results/6_classification_agreement.Rds"
saveAgreementFileCsv <- "../../data/BSD_results/6_classification_agreement.csv"

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
statisticsNames <- get_statistics_names(parNames)
designNames <- statisticsNames$design
statisticsNames <- statisticsNames[which(names(statisticsNames)!="design")]

#################################################
# Make test splits to cover all images
#################################################
nSegments <- length(unique(segmentStats$ImageName))
sampleSegments <- sample(unique(segmentStats$ImageName))
splitSize <- ceiling(nSegments/dataSplits)
testSegmentsList <- split(sampleSegments,
                          ceiling(seq_along(sampleSegments)/splitSize))

allDataTask <- make_task_BSD(segmentStats)

agreementDf <- NULL
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
  betterHOSPairs <- list()
  betterFASHOSPairs <- list()
  wrongPairsFAS <- list()
  wrongPairsHOS <- list()
  wrongPairsFASHOS <- list()
  correctPairsFAS <- list()
  correctPairsHOS <- list()
  correctPairsFASHOS <- list()
  for (s in c(1:2)) {
    trainData <- trainDataSplit[[s]]
    # fit to FAS
    outputFAS <- train_test_ridge(trainData, testData, statisticsNames$FAS,
                     balanceWeights=TRUE, subsetsPCA=NA)
    outputHOS <- train_test_ridge(trainData, testData, statisticsNames$HOS,
                     balanceWeights=TRUE, subsetsPCA=NA)
    outputFASHOS <- train_test_ridge(trainData, testData,
                                     unlist_names(statisticsNames[c(2,3)]),
                                     balanceWeights=TRUE, subsetsPCA=NA)

    ##### Get pairs where FAS failed but HOS or FAS+HOS succeeded
    wrongPairsFAS[[s]] <- which(outputFAS$correctPredictions==0)
    wrongPairsHOS[[s]] <- which(outputHOS$correctPredictions==0)
    wrongPairsFASHOS[[s]] <- which(outputFASHOS$correctPredictions==0)
    correctPairsFAS[[s]] <-  which(outputFAS$correctPredictions==1)
    correctPairsHOS[[s]] <-  which(outputHOS$correctPredictions==1)
    correctPairsFASHOS[[s]] <- which(outputFASHOS$correctPredictions==1)

    betterHOSPairs[[s]] <- wrongPairsFAS[[s]][which(wrongPairsFAS[[s]] %in%
                                               correctPairsHOS[[s]])]
    betterFASHOSPairs[[s]] <- wrongPairsFAS[[s]][which(wrongPairsFAS[[s]] %in%
                                                  correctPairsFASHOS[[s]])]
  }

  wrongFASAgree <- wrongPairsFAS[[1]][which(wrongPairsFAS[[1]] %in%
                                            wrongPairsFAS[[2]])]
  wrongHOSAgree <- wrongPairsHOS[[1]][which(wrongPairsHOS[[1]] %in%
                                            wrongPairsHOS[[2]])]
  wrongFASHOSAgree <- wrongPairsFASHOS[[1]][which(wrongPairsFASHOS[[1]] %in%
                                            wrongPairsFASHOS[[2]])]
  correctFASAgree <- correctPairsFAS[[1]][which(correctPairsFAS[[1]] %in%
                                            correctPairsFAS[[2]])]
  correctHOSAgree <- correctPairsHOS[[1]][which(correctPairsHOS[[1]] %in%
                                            correctPairsHOS[[2]])]
  correctFASHOSAgree <- correctPairsFASHOS[[1]][which(correctPairsFASHOS[[1]] %in%
                                            correctPairsFASHOS[[2]])]
  betterHOSAgree <- betterHOSPairs[[1]][which(betterHOSPairs[[1]] %in%
                                                 betterHOSPairs[[2]])]
  betterFASHOSAgree <- betterFASHOSPairs[[1]][which(betterFASHOSPairs[[1]] %in%
                                                       betterFASHOSPairs[[2]])]
  testData$betterHOS <- 0 
  testData$betterHOS[betterHOSAgree] <- 1
  testData$betterFASHOS <- 0 
  testData$betterFASHOS[betterFASHOSAgree] <- 1
  testData$wrongFAS <- 0 
  testData$wrongFAS[wrongFASAgree] <- 1
  testData$wrongHOS <- 0 
  testData$wrongHOS[wrongHOSAgree] <- 1
  testData$wrongFASHOS <- 0 
  testData$wrongFASHOS[wrongFASHOSAgree] <- 1
  testData$correctFAS <- 0 
  testData$correctFAS[correctFASAgree] <- 1
  testData$correctHOS <- 0 
  testData$correctHOS[correctHOSAgree] <- 1
  testData$correctFASHOS <- 0 
  testData$correctFASHOS[correctFASHOSAgree] <- 1

  agreementDf <- rbind(agreementDf, testData)
  # save dataframe with segment pairs where HOS is better
  saveRDS(agreementDf, saveAgreementFile)
  progressStr <- paste("Data split: ", r, "/", dataSplits, sep="")
  print(progressStr)
}
write.csv(agreementDf, saveAgreementFileCsv, row.names=FALSE)

