library(tidyr)
library(dplyr)
source("./analysis_functions.R")
set.seed(2691)

#############################
#############################
#############################
# Get the dimensionality of PCA for different
# statistics
#############################
#############################
#############################
pcaVar <- 0.95

#############################
# load data
#############################
dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
segmentStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.) %>%
  remove_constant_stats(.)

# extract names
parNames <- names(segmentStats)
statisticsNames <- get_statistics_names(parNames)
designNames <- statisticsNames$design
statisticsNames <- statisticsNames[which(names(statisticsNames)!="design")]

##############################################
# Do PCA for coarser subsets
##############################################
pairsStats <- make_task_BSD(segmentStats)

pixel <- c(0,1)
FAS <- c(0,1)
HOS <- c(0,1)
statsTypes <- c("pixel", "FAS", "HOS")
noSubsetMat <- expand.grid(pixel, FAS, HOS) %>%
  dplyr::mutate(., PCA_dim=NA, original_dim=NA) %>%
  dplyr::rename(., pixel=Var1, FAS=Var2, HOS = Var3) %>%
  dplyr::filter(., !(pixel==0 & FAS==0 & HOS==0))
resultsDf <- NULL

# No PCA subsetting
for (m in c(1:nrow(noSubsetMat))) {
  statsInd <- which(c(noSubsetMat[m,c("pixel", "FAS", "HOS")])==1)
  trialTypes <- statsTypes[statsInd]
  trialStatsList <- statisticsNames[trialTypes]
  trialStatsVec <- unlist_names(trialStatsList) 
  preparedData <- prepare_data_fit_test(trainData=pairsStats,
                        testData=pairsStats,
                        statsToUse=trialStatsVec,
                        balanceWeights=FALSE,
                        subsetsPCA=list(all=trialStatsVec),
                        labelColumn="same",
                        varianceRetained=pcaVar)
  pcaStats <- preparedData$trainStats
  noSubsetMat$PCA_dim[m] <- ncol(pcaStats)
  noSubsetMat$original_dim[m] <- length(trialStatsVec)
}

##############################################
# Do PCA for finer subsets
##############################################
statisticsNamesFiner <- get_statistics_names(parNames, subsetHOS=TRUE)
statisticsNamesFiner <- statisticsNamesFiner[which(names(statisticsNamesFiner)!="design")]
statisticsNamesFiner <- c(list(pixel=statisticsNamesFiner$pixel),
                          list(FAS=statisticsNamesFiner$FAS),
                          statisticsNamesFiner$HOS)

pixel <- c(0,1)
FAS <- c(0,1)
acm <- c(0,1)
cmc <- c(0,1)
pmc <- c(0,1)
prc <- c(0,1)
statsTypes <- c("pixel", "FAS", "acm", "cmc", "pmc", "prc")

subsetMat <- expand.grid(pixel, FAS, acm, cmc, pmc, prc) %>%
  dplyr::mutate(., PCA_dim=NA, original_dim=NA) %>%
  dplyr::rename(., pixel=Var1, FAS=Var2, acm=Var3, cmc=Var4,
                pmc=Var5, prc=Var6) %>%
  dplyr::filter(., (pixel+FAS+acm+cmc+pmc+prc)!=0)

for (m in c(1:nrow(subsetMat))) {
  statsInd <- which(c(subsetMat[m,c("pixel", "FAS", "acm", "cmc", "pmc",
                                        "prc")])==1)
  trialTypes <- statsTypes[statsInd]
  trialStatsList <- statisticsNamesFiner[trialTypes]
  trialStatsVec <- unlist_names(trialStatsList) 
  preparedData <- prepare_data_fit_test(trainData=pairsStats,
                        testData=pairsStats,
                        statsToUse=statsVec,
                        balanceWeights=FALSE,
                        subsetsPCA=list(all=trialStatsVec),
                        labelColumn="same",
                        varianceRetained=0.95)
  pcaStats <- preparedData$trainStats
  subsetMat$PCA_dim[m] <- ncol(pcaStats)
  subsetMat$original_dim[m] <- length(trialStatsVec)
}

singleSubset <- dplyr::filter(subsetMat, pixel+FAS+acm+cmc+pmc+prc==1)


#############################
#############################
#############################
# Run linear regression fixin number
# of PC for HOS
#############################
#############################
#############################
set.seed(2691)
dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
saveResults <- "../../data/BSD_results/3_BSD_results_PCA_60_matchedN.RDS"
repExp <- 20
pcaVar <- 0.60

#############################
# load data
#############################
segmentStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.) %>%
  remove_constant_stats(.)

parNames <- names(segmentStats)
statisticsNames <- get_statistics_names(parNames)
designNames <- statisticsNames$design
statisticsNames <- statisticsNames[which(names(statisticsNames)!="design")]

#############################
#generate template of design matrix for one repetition
#############################
pixel <- c(0,1)
FAS <- c(0,1)
HOS <- c(1)
statsTypes <- c("pixel", "FAS", "HOS")
designMatrixTemp <- expand.grid(pixel, FAS, HOS) %>%
  dplyr::mutate(., rep=NA, performance=NA) %>%
  dplyr::rename(., pixel=Var1, FAS=Var2, HOS=Var3) %>%
  dplyr::filter(., !(pixel==0 & FAS==0 & HOS==0))
resultsDf <- NULL

#############################
#Put together the pairs of patches
#############################
allDataTask <- make_task_BSD(segmentStats)

for (r in 1:repExp) {
  # split data into train and test set
  nSegments <- length(unique(segmentStats$ImageName))
  sampleSegments <- sample(unique(segmentStats$ImageName))
  trainSegments <- sampleSegments[1:floor(nSegments*4/5)]
  testSegments <- sampleSegments[(floor(nSegments*4/5)+1):nSegments]

  trainData <- dplyr::filter(allDataTask, ImageName %in% trainSegments) %>%
    droplevels(.)
  testData <- dplyr::filter(allDataTask, ImageName %in% testSegments) %>%
    droplevels(.)
  copyTemplate <- designMatrixTemp %>%
    dplyr::mutate(., rep=r)

  for (m in c(1:nrow(copyTemplate))) {
    statsInd <- which(c(copyTemplate[m,c("pixel", "FAS", "HOS")])==1)
    trialTypes <- statsTypes[statsInd]
    trialStatsList <- statisticsNames[trialTypes]
    trialStatsVec <- unlist_names(statisticsNames) 

    # Do PCA and keep only fixed number of HOS PC
    preparedData <- prepare_data_fit_test(trainData=trainData, testData=testData,
                          statsToUse=trialStatsVec,
                          balanceWeights=TRUE,
                          subsetsPCA=statisticsNames,
                          labelColumn="same",
                          varianceRetained=pcaVar)
    FAS_PC <- grep("FAS", names(preparedData$trainStats))
    HOS_PC <- grep("HOS", names(preparedData$trainStats))
    HOS_PC_rem <- HOS_PC[c((length(FAS_PC)+1):length(HOS_PC))]
    preparedData$trainStats <- preparedData$trainStats[,-HOS_PC_rem]
    preparedData$testStats <- preparedData$testStats[,-HOS_PC_rem]
    retainInd <- NULL
    trialNames <- names(statsInd)
    for (l in c(1:length(trialNames))) {
      retainInd <- c(retainInd, grep(trialNames[l],
                                           names(preparedData$trainStats)))
    }
    preparedData$trainStats <- preparedData$trainStats[,retainInd]
    preparedData$testStats <- preparedData$testStats[,retainInd]
    modelFit <- cv.glmnet(x=as.matrix(preparedData$trainStats),
                         y=preparedData$trainLabel,
                         weights=preparedData$weights,
                         family="binomial",
                         alpha=0, # alpha = 0 is ridge regression
                         standardize=FALSE,
                         type.measure="class")

    modelPredictions <- predict(modelFit, as.matrix(preparedData$testStats),
                                   type="class")
    correctPredictions <- as.integer(modelPredictions == preparedData$testLabel)
    predictionOutcome <- mean(correctPredictions)

    copyTemplate$performance[m] <- predictionOutcome
    print(paste("Rep: ", r,"/", repExp, "     Row: ", m, "/",
                nrow(copyTemplate), sep=""))
  }
  resultsDf <- rbind(resultsDf, copyTemplate)
  saveRDS(resultsDf, saveResults)
}


