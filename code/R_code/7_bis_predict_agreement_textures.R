library(tidyr)
library(dplyr)
library(glmnet)
library(caret)
library(vcd)
source("./analysis_functions.R")
set.seed(2691)

agreementFile <- "../../data/texture_results/6_classification_agreement_textures.Rds"
saveConfusionMat <- "../../data/texture_results/7_confusion_matrices_texture.RDS"
saveTextureOutcome <- "../../data/BSD_results/7_segment_outcomes.RDS"
dataSplits <- 10

#############################
# load data
#############################
agreementDf <- readRDS(agreementFile) %>%
  as_tibble(.) %>%
  remove_constant_stats(.)

#############################
# get the names of the different stats to use
#############################
parNames <- names(agreementDf)
statisticsNames <- get_statistics_names(parNames)
designNames <- statisticsNames$design
statisticsNames <- statisticsNames[which(names(statisticsNames)!="design")]

################
# Do linear regression to see whether the pairs with better HOS are clearly
# separated
################

# Make test splits to cover all images
nTextures <- length(unique(agreementDf$texture1))
sampleTextures <- sample(unique(agreementDf$texture1))
splitSize <- ceiling(nTextures/dataSplits)
testTexturesList <- split(sampleTextures,
                          ceiling(seq_along(sampleTextures)/splitSize))

confusionMatHOS <- matrix(0,2,2)
confusionMatFASHOS <- matrix(0,2,2)
texturesOutcome <- NULL

for (r in c(1:length(testTexturesList))) {
  testTextures <- testTexturesList[[r]]
  trainTextures <- sampleTextures[which(!sampleTextures %in% testTextures)]

  trainData <- dplyr::filter(agreementDf, texture1 %in% trainTextures) %>%
    droplevels(.)
  testData <- dplyr::filter(agreementDf, texture1 %in% testTextures) %>%
    droplevels(.)

  statsToUse <- unlist_names(statisticsNames[c("FAS", "HOS")])

  predictBetterHOS <- train_test_ridge(trainData, testData,
                                  statsToUse=statsToUse,
                                  balanceWeights=TRUE, subsetsPCA=NA,
                                  labelColumn="betterHOS")

  predictBetterFASHOS <- train_test_ridge(trainData, testData,
                                  statsToUse=statsToUse,
                                  balanceWeights=TRUE, subsetsPCA=NA,
                                  labelColumn="betterFASHOS")

  # add results to global confusion matrix
  confusionMatHOS <- confusionMatHOS + predictBetterHOS$confusionMatrix$table
  confusionMatFASHOS <- confusionMatFASHOS + predictBetterFASHOS$confusionMatrix$table

  # save in a Df the status of each pair (false positive, false negative, etc)
  predictOutcomeHOS <- confusion_indices(predictBetterHOS$predictions,
                                    as.character(testData$betterHOS))
  predictOutcomeFASHOS <- confusion_indices(predictBetterFASHOS$predictions,
                                    as.character(testData$betterFASHOS))

  shorterTest <- dplyr::select(testData, designNames) %>%
    dplyr::mutate(., predictionHOS=predictOutcomeHOS,
                  predictionFASHOS=predictOutcomeFASHOS)
  texturesOutcome <- rbind(texturesOutcome, shorterTest)   

  
}

confusionMatrices <- list(confusionHOS=confusionMatHOS,
                          confusionFASHOS=confusionMatFASHOS)

saveRDS(confusionMatrices, saveConfusionMat)
saveRDS(texturesOutcome, saveTextureOutcome)

