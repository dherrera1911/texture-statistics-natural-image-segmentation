library(tidyr)
library(dplyr)
library(glmnet)
library(matrixStats)
source("./analysis_functions.R")
set.seed(2691)

dataFile <- "../../data/texture_stats/texture_stats_statsNorm.csv"
saveResults <- "../../data/texture_results/texture_dnn.Rds"

repExp <- 20
nRep <- 5
layerUnits <- c(30, 10)
regularizationWeight <- 0.002
epochs <- 200


#############################
# load data
#############################
textureStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.)

#############################
# get the names of the different stats to use
#############################
parNames <- names(textureStats)
statisticsNames <- get_statistics_names(parNames)
designNames <- statisticsNames$design
statisticsNames <- statisticsNames[which(names(statisticsNames)!="design")]

#############################
#generate template of design matrix for one repetition
#############################
pixel <- c(0,1)
FAS <- c(0,1)
HOS <- c(0,1)
statsTypes <- c("pixel", "FAS", "HOS")
designMatrixTemp <- expand.grid(pixel, FAS, HOS) %>%
  dplyr::mutate(., rep = NA, performance = NA) %>%
  dplyr::rename(., pixel = Var1, FAS = Var2, HOS = Var3) %>%
  dplyr::filter(., !(pixel==0 & FAS==0 & HOS==0))
resultsDf <- NULL


for (r in 1:repExp) {
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

  copyTemplate <- designMatrixTemp %>%
    dplyr::mutate(., rep = r)

  for (m in c(1:length(nrow(copyTemplate)))) {
    statsInd <- which(c(copyTemplate[m,c("pixel", "FAS", "HOS")])==1)
    trialTypes <- statsTypes[statsInd]
    trialStatsList <- statisticsNames[trialTypes]
    trialStatsVec <- unlist_names(trialStatsList) 
    modelOutcome <- train_test_dnn(trainData=trainData, testData=testData,
                     statsToUse=trialStatsVec, balanceWeights=TRUE,
                     subsetsPCA=list(all=trialStatsVec),
                     layerUnits=layerUnits,
                     regularizationWeight=regularizationWeight,
                     epochs=epochs)
    copyTemplate$performance[m] <- modelOutcome$accuracy
  }
  resultsDf <- rbind(resultsDf, copyTemplate)
}

saveRDS(resultsDf, saveResults)



