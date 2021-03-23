library(tidyr)
library(dplyr)
library(glmnet)
source("./analysis_functions.R")
set.seed(2691)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_pixNorm.csv"
dataFile <- "../../data/texture_stats/texture_stats_statsNorm.csv"
savePerformanceResults <- "../../data/texture_results/8_subset_stats_performance_texture.RDS"
#saveTexturePredictions <- "../../data/BSD_results/8_bsd_predictions_subsets.RDS"
#dataSplits <- 10
repExp <- 20
nRep <- 5

#############################
# load data
#############################
textureStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.) %>%
  remove_constant_stats(.)

#############################
# get the names of the different stats to use
#############################
parNames <- names(textureStats)
statisticsNames <- get_statistics_names(parNames, subsetHOS=TRUE)
designNames <- statisticsNames$design
statisticsNames <- statisticsNames[which(names(statisticsNames)!="design")]

#############################
#generate template of design matrix for one repetition
#############################
HOStypes <- names(statisticsNames$HOS)
FAS <- c(0,1)
acm <- c(0,1)
cmc <- c(0,1)
pmc <- c(0,1)
prc <- c(0,1)

designMatrixTemp <- expand.grid(FAS, acm, cmc, pmc, prc) %>%
  dplyr::mutate(., rep=NA, performance=NA) %>%
  dplyr::rename(., FAS=Var1, acm=Var2, cmc=Var3,
                pmc=Var4, prc=Var5) %>%
  dplyr::filter(., (FAS+acm+cmc+pmc+prc)!=0)

modelPerformance <- NULL
texturePredictions <- NULL

#############################
#Fit the models
#############################
#for (r in 1:dataSplits) {
for (r in 1:repExp) {
  # sample train and test data examples
  nTextures <- length(unique(textureStats$texture))
  sampleTextures <- sample(unique(textureStats$texture))
  trainTextures <- sampleTextures[1:round(nTextures/2)]
  testTextures <- sampleTextures[(1+round(nTextures/2)):nTextures]

  trainData <- dplyr::filter(textureStats, texture %in% trainTextures) %>%
    droplevels(.) %>%
    make_task_textures(., nRep)
  testData <- dplyr::filter(textureStats, texture %in% testTextures) %>%
    droplevels(.) %>%
    make_task_textures(., nRep)

  copyTemplate <- designMatrixTemp %>%
    dplyr::mutate(., rep = r)

  parNames <- names(trainData)
  statisticsNames <- get_statistics_names(parNames, subsetHOS=TRUE)
  designNames <- statisticsNames$design

  textureResultsDf <- dplyr::select(testData, same, all_of(designNames))

  for (m in c(1:nrow(copyTemplate))) {
    statsInd <- which(c(copyTemplate[m,c("acm", "cmc", "pmc", "prc")])==1)
    trialHOS <- HOStypes[statsInd]
    if (length(trialHOS)>0) {
      trialHOSList <- statisticsNames$HOS[trialHOS]
      trialStatsVec <- unlist_names(trialHOSList) 
      modelName <- paste(trialHOS, collapse="_")
    } else {
      trialStatsVec <- NULL
      modelName <- NULL
    }
    if (copyTemplate$FAS[m] == 1) {
      trialStatsVec <- c(statisticsNames$FAS, trialStatsVec)
      if (length(trialHOS) > 0) {
        modelName <- paste("FAS_", modelName, sep="")
      } else {
        modelName <- "FAS"
      }
    }
    print(paste("Row: ", m, "/", nrow(copyTemplate), sep=""))
    modelOutcome <- train_test_ridge(trainData=trainData, testData=testData,
                     statsToUse=trialStatsVec, balanceWeights=TRUE)
    copyTemplate$performance[m] <- modelOutcome$accuracy
#    textureResultsDf[[modelName]] <- modelOutcome$predictions
  }
  modelPerformance <- rbind(modelPerformance, copyTemplate)
#  texturePredictions <- rbind(texturePredictions, textureResultsDf)
  saveRDS(modelPerformance, savePerformanceResults)
#  saveRDS(texturePredictions, saveTexturePredictions)
}

