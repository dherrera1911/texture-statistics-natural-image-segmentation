library(tidyr)
library(dplyr)
library(ggplot2)
library(matrixStats)
source("./analysis_functions.R")
set.seed(2691)
#library(Hmisc)

plottingDir <- "../../data/plots/"

##################################
##################################
##### 1 Plot N of parameters #####
##################################
##################################

statsFile <- "../../data/texture_stats/texture_stats_statsNorm.csv"
statsDf <- read.csv(statsFile, sep = ",") %>%
  as_tibble(.) %>%
  remove_constant_stats(.)

# get the names of the different stats to use
parNames <- names(statsDf)
statisticsNames <- get_statistics_names(parNames)
designNames <- statisticsNames$design

nPix <- length(statisticsNames$pixel)
nFAS <- length(statisticsNames$FAS)
nHOS <- length(statisticsNames$HOS)

statFactor <- factor(c("Pix", "FAS", "HOS"), levels = c("Pix", "FAS", "HOS"))
nParDf <- data.frame(Statistics = statFactor,
                     nPar = c(nPix, nFAS, nHOS))

parPlot <- dplyr::filter(nParDf, Statistics %in% c("FAS", "HOS")) %>%
  ggplot(data=., aes(x = Statistics, y = nPar, fill = Statistics)) +
  geom_bar(stat="identity", width = 0.5) +
  scale_fill_manual(name = "Statistic groups",
                     #labels = c("Pixel", "Spectral", "HOS"),
                     #  values = c("#000000", "#2691d4", "#ba2229"),
                     labels = c("Spectral", "HOS"),
                       values = c("#2691d4", "#ba2229"),
                     guide = "none") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid")) +
  ylab("Number of parameters") +
  scale_x_discrete(name = "Statistics",
                   #labels = c("Pixel", "FAS", "HOS"))
                   labels = c("Spectral", "HOS"))

ggsave(paste(plottingDir, "nParameters.png", sep=""), parPlot,
       width = 5, height = 5, units = "cm") 

#############################
#############################
##### 2 Texture results #####
#############################
#############################

assign_stat_name <- function(pixel, FAS, HOS) {
  name <- NULL
  for (i in c(1:length(pixel))) {
    pixel_i <- pixel[i]
    FAS_i <- FAS[i]
    HOS_i <- HOS[i]
    if (pixel_i & !FAS_i & !HOS_i) {
      name[i] <- "Pix"
    } else if (pixel_i & FAS_i & !HOS_i) {
      name[i] <- "Pix_FAS"
    } else if (pixel_i & !FAS_i & HOS_i) {
      name[i] <- "Pix_HOS"
    } else if (pixel_i & FAS_i & HOS_i) {
      name[i] <- "Pix_FAS_HOS"
    } else if (!pixel_i & FAS_i & HOS_i) {
      name[i] <- "FAS_HOS"
    } else if (!pixel_i & FAS_i & !HOS_i) {
      name[i] <- "FAS"
    } else if (!pixel_i & !FAS_i & HOS_i) {
      name[i] <- "HOS"
    }
  }
  return(name)
}

assign_stat_name_plot <- function(pixel, FAS, HOS) {
  name <- NULL
  for (i in c(1:length(pixel))) {
    pixel_i <- pixel[i]
    FAS_i <- FAS[i]
    HOS_i <- HOS[i]
    if (pixel_i & !FAS_i & !HOS_i) {
      name[i] <- "Pix"
    } else if (pixel_i & FAS_i & !HOS_i) {
      name[i] <- "Pix-Spectral"
    } else if (pixel_i & !FAS_i & HOS_i) {
      name[i] <- "Pix-HOS"
    } else if (pixel_i & FAS_i & HOS_i) {
      name[i] <- "Pix-Spectral-HOS"
    } else if (!pixel_i & FAS_i & HOS_i) {
      name[i] <- "Spectral-HOS"
    } else if (!pixel_i & FAS_i & !HOS_i) {
      name[i] <- "Spectral"
    } else if (!pixel_i & !FAS_i & HOS_i) {
      name[i] <- "HOS"
    }
  }
  return(name)
}

textureResultsFile <- "../../data/texture_results/texture_results.Rds"
textureRes <- readRDS(textureResultsFile) %>%
  dplyr::mutate(., statsName=assign_stat_name_plot(pixel, FAS, HOS),
                error=1-performance)
textureRes$statsName <- factor(textureRes$statsName,
                               levels=c("Pix", "Pix-Spectral",
                                        "Pix-HOS", "Pix-Spectral-HOS",
                                        "Spectral", "HOS", "Spectral-HOS"))

texturePerformanceSummary <- group_by(textureRes, statsName, pixel) %>%
  summarize(., meanError=mean(error), sdError=sd(error))

#### plot results when using pixel data ####
texturePlot <- dplyr::filter(textureRes, pixel==1) %>%
  ggplot(aes(x=statsName, y=error*100, color=statsName)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, shape = 1) +
#  stat_summary(fun = mean_cl_normal, geom = "errorbar") +
  #geom_boxplot() +
  scale_color_manual(name = "Statistic groups",
                     labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS"),
                       values = c("#000000", "#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  stat_summary(fun = "mean") +
  ylim(0, 30) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid")) +
  ylab("Error rate (%)") +
  scale_x_discrete(name = "Parameters") +
  NULL

texturePlotName <- paste(plottingDir, "texturePlot.png", sep="")
ggsave(texturePlotName, texturePlot, width = 10, height = 7, units = "cm") 

#### plot results without pixel data ####
texturePlot2 <- dplyr::filter(textureRes, pixel==0) %>%
  ggplot(aes(x=statsName, y=error*100, color=statsName)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, shape = 1) +
#  stat_summary(fun = mean_cl_normal, geom = "errorbar") +
  #geom_boxplot() +
  scale_color_manual(name = "Statistic groups",
                     labels = c("Spectral", "HOS", "Spectral-HOS"),
                       values = c("#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  stat_summary(fun = "mean") +
  ylim(0, 30) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid")) +
  ylab("Error rate (%)") +
  scale_x_discrete(name = "Parameters") +
  NULL

texturePlotName <- paste(plottingDir, "texturePlot_noPix.png", sep="")
ggsave(texturePlotName, texturePlot2, width = 10, height = 7, units = "cm") 


##########################
##### 3 Texture DNN ######
##########################

#textureResultsDNNFile <- "../../data/texture_results/texture_dnn.Rds"
#
#textureDNNRes <- readRDS(textureResultsDNNFile) %>%
#  pivot_longer(., c("Pix", "FA", "HOS", "FA_HOS"), names_to = "Model_params",
#               values_to = "Performance") %>%
#  dplyr::filter(., !is.na(Performance))
#
## add clearer the stats
#colFAS <- rep(0, nrow(textureDNNRes))
#colFAS[grep("FA", textureDNNRes$Model_params)] <- 1
#colHOS <- rep(0, nrow(textureDNNRes))
#colHOS[grep("*HOS", textureDNNRes$Model_params)] <- 1
#
#textureDNNRes$FAS <- colFAS
#textureDNNRes$HOS <- colHOS
#
#pixRes <- dplyr::filter(textureDNNRes, usePix == 1) %>%
#  dplyr::mutate(., Model_params =
#                c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS")[usePix + FAS + 2*HOS])
#
#modelOrder <- c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS", "FA", "HOS", "FA_HOS")
#
#textureDNNRes <- dplyr::filter(textureDNNRes, usePix == 0) %>%
#  rbind(., pixRes) %>%
#  dplyr::mutate(., Error = 1 - Performance)
#
#textureDNNRes$Model_params <- factor(textureDNNRes$Model_params, levels = modelOrder)
#
#texturePerformanceDNN <- textureDNNRes %>%
#  dplyr::select(., Model_params, Error) %>%
#  group_by(., Model_params) %>%
#  dplyr::summarize(., mean_error = mean(Error), error_sd = sd(Error))
#
#
##### plot results when using pixel data ####
#textureDNNPlot <- dplyr::filter(textureDNNRes, usePix == 1) %>%
#  ggplot(aes(x = Model_params, y = Error  * 100, color = Model_params)) +
#  geom_jitter(width = 0.2, height = 0, alpha = 0.5, shape = 1) +
#  scale_color_manual(name = "Statistic groups",
#                     labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS"),
#                       values = c("#000000", "#2691d4", "#b400e5", "#ba2229"),
#                     guide = "none") +
#  stat_summary(fun = "mean") +
#  ylim(0, 30) +
#  theme_bw() +
#  theme(panel.border = element_blank(),
#        axis.line.x = element_line(size=0.5, linetype="solid"),
#        axis.line.y = element_line(size=0.5, linetype="solid")) +
#  ylab("Error rate (%)") +
#  scale_x_discrete(name = "Parameters",
#                   labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS")) +
#  NULL
#
#textureDNNPlotName <- paste(plottingDir, "textureDNNPlot.png", sep="")
#ggsave(textureDNNPlotName, textureDNNPlot, width = 10, height = 7, units = "cm") 
#
#

#############################
#############################
##### 4 BSD results #########
#############################
#############################

bsdFile <- "../../data/BSD_results/3_BSD_results.RDS"
bsdRes <- readRDS(bsdFile) %>%
  dplyr::mutate(., statsName=assign_stat_name_plot(pixel, FAS, HOS),
                error=1-performance)
bsdRes$statsName <- factor(bsdRes$statsName,
                               levels=c("Pix", "Pix-Spectral",
                                        "Pix-HOS", "Pix-Spectral-HOS",
                                        "Spectral", "HOS", "Spectral-HOS"))

bsdPerformanceSummary <- group_by(bsdRes, statsName, pixel) %>%
  summarize(., meanError=mean(error), sdError=sd(error))

## plot results when using pixel data ####
bsdPlot <- dplyr::filter(bsdRes, pixel==1) %>%
  droplevels(.) %>%
  ggplot(aes(x = statsName, y = error  * 100, color = statsName)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, shape = 1) +
  scale_color_manual(name = "Statistic groups",
                     labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS"),
                       values = c("#000000", "#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  stat_summary(fun = "mean") +
  ylim(0, 30) +
  theme_bw() +
  ylab("error rate (%)") +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid")) +
  scale_x_discrete(name = "Parameters")

## plot results when using pixel data ####
bsdPlot2 <- dplyr::filter(bsdRes, pixel==0) %>%
  droplevels(.) %>%
  ggplot(aes(x = statsName, y = error  * 100, color = statsName)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.5, shape = 1) +
  scale_color_manual(name = "Statistic groups",
                     labels = c("Spectral", "HOS", "Spectral-HOS"),
                       values = c("#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  stat_summary(fun = "mean") +
  ylim(0, 30) +
  theme_bw() +
  ylab("error rate (%)") +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid")) +
  scale_x_discrete(name = "Parameters")

BSDPlotName <- paste(plottingDir, "BSDPlot_noPix.png", sep="")
ggsave(BSDPlotName, bsdPlot2, width = 10, height = 7, units = "cm") 

#########################
#########################
##### 5 BSD DNN #########
#########################
#########################

bsdDNNFile <- "../../data/BSD_results/4_BSD_dnn_single.RDS"
bsdDNNRes <- readRDS(bsdDNNFile) %>%
  dplyr::mutate(., statsName=assign_stat_name_plot(pixel, FAS, HOS),
                error=1-performance)
bsdDNNRes$statsName <- factor(bsdDNNRes$statsName,
                               levels=c("Pix", "Pix-Spectral",
                                        "Pix-HOS", "Pix-Spectral-HOS",
                                        "Spectral", "HOS", "Spectral-HOS"))

bsdDNNPerformanceSummary <- group_by(bsdDNNRes, statsName, pixel,
                                     architecture, regularization, epochs) %>%
  summarize(., meanError=mean(error), sdError=sd(error), count=n())

bsdDNNPlot <- dplyr::filter(bsdDNNRes, pixel==1) %>%
  droplevels(.) %>%
  ggplot(aes(x=statsName, y=error*100, color=statsName)) +
  geom_jitter(width=0.2, height=0, alpha=0.5, shape=1) +
  scale_color_manual(name = "Statistic groups",
                       values = c("#000000", "#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  stat_summary(fun = "mean") +
  ylim(0, 30) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid")) +
  ylab("error rate (%)") +
  scale_x_discrete(name = "Parameters")

bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot.png", sep="")
ggsave(bsdDNNPlotName, bsdDNNPlot, width = 10, height = 7, units = "cm") 

bsdDNNPlot2 <- dplyr::filter(bsdDNNRes, pixel==0) %>%
  droplevels(.) %>%
  ggplot(aes(x=statsName, y=error*100, color=statsName)) +
  geom_jitter(width=0.2, height=0, alpha=0.5, shape=1) +
  scale_color_manual(name = "Statistic groups",
                       values = c("#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  stat_summary(fun="mean") +
  ylim(0, 30) +
  theme_bw() +
  theme(panel.border=element_blank(),
        axis.line.x=element_line(size=0.5, linetype="solid"),
        axis.line.y=element_line(size=0.5, linetype="solid")) +
  ylab("error rate (%)") +
  scale_x_discrete(name = "Parameters")

bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot_noPix.png", sep="")
ggsave(bsdDNNPlotName, bsdDNNPlot2, width = 10, height = 7, units = "cm") 

##################################
##################################
##### 5.2 BSD DNN params #########
##################################
##################################

bsdDNNFile <- "../../data/BSD_results/4_BSD_dnn.RDS"
bsdDNNRes <- readRDS(bsdDNNFile) %>%
  dplyr::mutate(., statsName=assign_stat_name_plot(pixel, FAS, HOS),
                error=1-performance)
bsdDNNRes$statsName <- factor(bsdDNNRes$statsName,
                               levels=c("Pix", "Pix-Spectral",
                                        "Pix-HOS", "Pix-Spectral-HOS",
                                        "Spectral", "HOS", "Spectral-HOS"))

bsdDNNPerformanceSummary <- group_by(bsdDNNRes, statsName, pixel,
                                     architecture, regularization, epochs) %>%
  summarize(., meanError=mean(error), sdError=sd(error), count=n())

# compare performance of NN to linear model
linearSpectralHOS <- dplyr::filter(bsdPerformanceSummary,
                                   statsName=="Pix-Spectral-HOS")[["meanError"]]

performanceDiff <- dplyr::filter(bsdDNNPerformanceSummary,
                                 statsName=="Pix-Spectral-HOS") %>%
  dplyr::mutate(., errorDiff=linearSpectralHOS-meanError)

bestFASHOS <- performanceDiff[which.min(performanceDiff$errorDiff),]


#### compare NN regimes
allEpochs <- unique(bsdDNNRes$epochs)
allReguls <- unique(bsdDNNRes$regularization)
bsdDNNParPlot <- list()
bsdDNNParPlot2 <- list()
for (ep in c(1:length(allEpochs))) {
  bsdDNNParPlot[[ep]] <- list()
  bsdDNNParPlot2[[ep]] <- list()
  for (reg in c(1:length(allReguls))) {
    bsdDNNParPlot[[ep]][[reg]] <- dplyr::filter(bsdDNNRes, pixel==1,
                                       epochs==allEpochs[ep] &
                                         regularization==allReguls[reg]) %>%
      droplevels(.) %>%
      ggplot(aes(x=statsName, y=error*100, color=statsName)) +
      geom_jitter(width=0.2, height=0, alpha=0.5, shape=1) +
      scale_color_manual(name = "Statistic groups",
                           values = c("#000000", "#2691d4", "#b400e5", "#ba2229"),
                         guide = "none") +
      stat_summary(fun = "mean") +
      ylim(0, 30) +
      theme_bw() +
      facet_wrap(facets=c("architecture", "regularization")) +
      theme(panel.border = element_blank(),
            axis.line.x = element_line(size=0.5, linetype="solid"),
            axis.line.y = element_line(size=0.5, linetype="solid")) +
      ylab("error rate (%)") +
      scale_x_discrete(name = "Parameters",
                       labels = c("Pix", "Pix-Spectral", "Pix-HOS",
                                  "Pix-Spectral-HOS"))
    bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot_param_", 
                            allEpochs[ep], "_", allReguls[reg],
                            ".png", sep="")
    ggsave(bsdDNNPlotName, bsdDNNParPlot[[ep]][[reg]],
           width=20, height=14, units = "cm") 
    bsdDNNParPlot2[[ep]][[reg]] <- dplyr::filter(bsdDNNRes, pixel==0,
                                       epochs==allEpochs[ep] &
                                         regularization==allReguls[reg]) %>%
      droplevels(.) %>%
      ggplot(aes(x=statsName, y=error*100, color=statsName)) +
      geom_jitter(width=0.2, height=0, alpha=0.5, shape=1) +
      scale_color_manual(name = "Statistic groups",
                           values = c("#2691d4", "#b400e5", "#ba2229"),
                         guide = "none") +
      stat_summary(fun = "mean") +
      ylim(0, 30) +
      theme_bw() +
      facet_wrap(facets=c("architecture", "regularization")) +
      theme(panel.border = element_blank(),
            axis.line.x = element_line(size=0.5, linetype="solid"),
            axis.line.y = element_line(size=0.5, linetype="solid")) +
      ylab("error rate (%)") +
      scale_x_discrete(name = "Parameters")
    bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot_param_", 
                            allEpochs[ep], "_", allReguls[reg],
                            "_noPix.png", sep="")
    ggsave(bsdDNNPlotName, bsdDNNParPlot2[[ep]][[reg]],
           width=20, height=14, units = "cm") 
  }
}


##################################
##################################
##### 6 analyze agreement #########
##################################
##################################
agreementFile <- "../../data/BSD_results/6_classification_agreement.csv"
agreementDf <- read.csv(agreementFile, stringsAsFactors=FALSE)

agreementSel <- dplyr::select(agreementDf, same, betterHOS, betterFASHOS) %>%
  tidyr::pivot_longer(., cols=c("betterHOS", "betterFASHOS"),
                      names_to="comparison", values_to="result")

agreementSummary <- agreementSel %>%
  group_by(., same, comparison) %>% 
  summarize(., Better=mean(result), Not_better=mean(result==0)) %>%
  tidyr::pivot_longer(., cols=c("Better", "Not_better"),
                      names_to="outcome", values_to="proportion") %>%
  dplyr::mutate(., comparison=factor(comparison,
                                       levels=c("betterHOS", "betterFASHOS"),
                                       labels=c("HOS", "Spectral-HOS")),
                outcome=factor(outcome,
                               levels=c("Not_better", "Better")),
                same=factor(same, levels=c(0, 1),
                            labels=c("Different", "Same")))

agreementPlot <- ggplot(data=agreementSummary, aes(x=outcome, y=proportion*100,
                                                   fill=same)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~comparison) +
  scale_fill_manual(name="Pair class",
                     labels=c("Different", "Same"),
                       values=c("#b33018", "#14b74b")) +
  xlab("") +
  ylab("% relative to spectral") +
  theme_bw()

agreementPlotName <- paste(plottingDir, "agreement_bars.png", sep="")
ggsave(agreementPlotName,  agreementPlot, width=10, height=6, units = "cm") 

agreementFASHOSPlot <- dplyr::filter(agreementSummary, comparison=="Spectral-HOS") %>%
  ggplot(data=., aes(x=outcome, y=proportion*100, fill=same)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(name="Pair class",
                     labels=c("Different", "Same"),
                       values=c("#b33018", "#14b74b")) +
  xlab("") +
  ylab("Compared to spectral") +
  theme_bw()

agreementPlotName <- paste(plottingDir, "agreement_bars_FASHOS.png", sep="")
ggsave(agreementPlotName,  agreementFASHOSPlot, width=10, height=6, units = "cm") 


##################################
##### 7 predict agreement ########
##################################
agreementPredFile <- "../../data/BSD_results/7_confusion_matrices.RDS"
agreementPredDf <- readRDS(agreementPredFile)

predictionHOS <- caret::confusionMatrix(agreementPredDf$confusionHOS)
predictionFASHOS <- caret::confusionMatrix(agreementPredDf$confusionFASHOS)


######################################
##### 8 analyze stats subsets ########
######################################
assign_subset_name_plot <- function(FAS, acm, cmc, pmc, prc) {
  name <- NULL
  for (i in c(1:length(FAS))) {
    FAS_i <- FAS[i]
    acm_i <- acm[i]
    cmc_i <- cmc[i]
    pmc_i <- pmc[i]
    prc_i <- prc[i]
    if (FAS_i) {
      name[i] <- "Spectral"
    } else {
      name[i] <- ""
    }
    if (acm_i) {
      name[i] <- paste(name[i], "Space", sep="_")
    }
    if (cmc_i) {
      name[i] <- paste(name[i], "Orientation", sep="_")
    }
    if (pmc_i) {
      name[i] <- paste(name[i], "Scale", sep="_")
    }
    if (prc_i) {
      name[i] <- paste(name[i], "Phase", sep="_")
    }
    if (strtrim(name[i],1)=="_") {
      nameTemp <- strsplit(name[i], "^_")
      name[i] <- nameTemp[[1]][2]
    }
  }
  return(name)
}

subsetsFile <- "../../data/BSD_results/8_subset_stats_performance_FAS.RDS"
subsetsDf <- readRDS(subsetsFile) %>%
  dplyr::mutate(., statsName=assign_subset_name_plot(FAS, acm, cmc, pmc, prc),
                error=1-performance,
                numberSubsets=FAS+acm+cmc+pmc+prc)

subsetsSummary <- group_by(subsetsDf, statsName, FAS, acm, cmc, pmc, prc, numberSubsets) %>%
  summarize(., meanError=mean(error), sdError=sd(error)/sqrt(n())) %>%
  dplyr::mutate(., numberSubsets=factor(numberSubsets))
subsetsSummary <- subsetsSummary[order(subsetsSummary$numberSubsets),]

fullModelVals <- dplyr::filter(subsetsSummary, numberSubsets %in% c(1,5))$meanError

addedSubsetsDf <- dplyr::filter(subsetsSummary, numberSubsets %in% c(2))
addedSubsetsDf$statsName <- sub("Spectral_", "+ ", addedSubsetsDf$statsName )

subsetsPlot <- addedSubsetsDf %>%
  ggplot(data=., aes(x=statsName, y=meanError*100)) +
  geom_bar(stat="identity", width = 0.5) +
  geom_errorbar(aes(ymin=(meanError-sdError*2)*100,
                    ymax=(meanError+sdError*2)*100), width=0.1) +
  geom_hline(yintercept=fullModelVals[1]*100, color="#2691d4", linetype=2) +
  geom_hline(yintercept=fullModelVals[2]*100, color="#ba2229", linetype=2) +
  ylab("error rate (%)") +
  xlab("Statistics used") +
  ylim(c(0, 18)) +
  theme_bw()

subsetsPlotName <- paste(plottingDir, "subsetsPlot_added.png", sep="")
ggsave(subsetsPlotName, subsetsPlot, width = 10, height = 7, units = "cm") 

removedSubsetsDf <- dplyr::filter(subsetsSummary, numberSubsets %in% c(4))
statsNames <- c("Space", "Orientation", "Scale", "Phase")
indVec <- NULL
nameVec <- NA
for (sN in statsNames) {
  ind <- which(!(c(1:4) %in% grep(sN, removedSubsetsDf$statsName)))
  nameVec[ind] <- sN
}
removedSubsetsDf$statsName <- paste("-", nameVec)

subsetsPlot2 <- removedSubsetsDf %>%
  ggplot(data=., aes(x=statsName, y=meanError*100)) +
  geom_bar(stat="identity", width = 0.5) +
  geom_errorbar(aes(ymin=(meanError-sdError*2)*100,
                    ymax=(meanError+sdError*2)*100), width=0.1) +
  geom_hline(yintercept=fullModelVals[1]*100, color="#2691d4", linetype=2) +
  geom_hline(yintercept=fullModelVals[2]*100, color="#ba2229", linetype=2) +
  ylab("error rate (%)") +
  xlab("Statistics used") +
  ylim(c(0, 18)) +
  theme_bw()

subsetsPlotName2 <- paste(plottingDir, "subsetsPlot_removed.png", sep="")
ggsave(subsetsPlotName2, subsetsPlot2, width = 10, height = 7, units = "cm") 

######################################
##### 9 Experiment results  ##########
######################################
expResultsFile <- "../../data/experiment_stats/experimentResults.RDS"
expSegmentationFile <- "../../data/BSD_results/9_segmentation_predictions.RDS"
expAgreementFile <- "../../data/BSD_results/9_better_stats_predictions.RDS"

expResults <- readRDS(expResultsFile) %>%
  as_tibble(.)

expSeg <- readRDS(expSegmentationFile) %>%
  as_tibble(.)

summaryExpSeg <- group_by(expSeg, texture) %>%
  summarize(., mean_same_FAS=mean(same_FAS),
            mean_same_HOS=mean(same_HOS),
            mean_same_FASHOS=mean(same_FASHOS),
            act_diff_FAS=-mean(activation_FAS),
            act_diff_HOS=-mean(activation_HOS),
            act_diff_FASHOS=-mean(activation_FASHOS),
            act_diff_FAS_sd=sd(activation_FAS)/sqrt(n()),
            act_diff_HOS_sd=sd(activation_HOS)/sqrt(n()),
            act_diff_FASHOS_sd=sd(activation_FASHOS)/sqrt(n()))

expAgreement <- readRDS(expAgreementFile) %>%
  as_tibble(.)

summaryExpAgreement <- group_by(expAgreement, texture) %>%
  summarize(., mean_betterHOS=mean(betterHOS),
            mean_betterFASHOS=mean(betterFASHOS),
            act_better_HOS=mean(betterHOSActivation),
            act_better_FASHOS=mean(betterFASHOSActivation),
            act_better_HOS_sd=sd(betterHOSActivation)/sqrt(n()),
            act_better_FASHOS_sd=sd(betterFASHOSActivation)/sqrt(n()))

# put all the data together
expDataFull <- merge(expResults, summaryExpSeg, by="texture") %>%
  merge(., summaryExpAgreement, by="texture") %>%
  dplyr::mutate(., texture=factor(texture, labels=c("T1", "T2", "T3", "T4")),
                estimate=-estimate, conf.low=-conf.low, conf.high=-conf.high)

experimentHOSSegPlot <- ggplot(expDataFull, aes(x=estimate, y=act_diff_HOS,
                                                color=texture)) +
  geom_point() +
  geom_errorbar(aes(ymin=act_diff_HOS-2*act_diff_HOS_sd,
                    ymax=act_diff_HOS+2*act_diff_HOS_sd)) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high)) +
  xlab("Experimental HOS segmentation") +
  ylab("Activation of segmentation model") +
  theme_bw() +
  NULL
experimentPlotName <- paste(plottingDir, "experimentModel_HOS.png", sep="")
ggsave(experimentPlotName, experimentHOSSegPlot, width=12, height=9, units = "cm") 

experimentFASHOSSegPlot <- ggplot(expDataFull, aes(x=estimate, y=act_diff_FASHOS,
                                                color=texture)) +
  geom_point() +
  geom_errorbar(aes(ymin=act_diff_FASHOS-2*act_diff_FASHOS_sd,
                    ymax=act_diff_FASHOS+2*act_diff_FASHOS_sd)) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high)) +
  xlab("Experimental HOS segmentation") +
  ylab("Activation of model segmentation") +
  theme_bw() +
  NULL
experimentPlotName <- paste(plottingDir, "experimentModel_FASHOS.png", sep="")
ggsave(experimentPlotName, experimentFASHOSSegPlot, width=12, height=9, units = "cm") 

experimentHOSAgreementPlot <- ggplot(expDataFull, aes(x=estimate, y=act_better_HOS,
                                                color=texture)) +
  geom_point() +
  geom_errorbar(aes(ymin=act_better_HOS-2*act_better_HOS_sd,
                    ymax=act_better_HOS+2*act_better_HOS_sd)) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high)) +
  xlab("Experimental HOS segmentation") +
  ylab("Activation of HOS improvement classifier") +
  theme_bw()

experimentFASHOSAgreementPlot <- ggplot(expDataFull, aes(x=estimate, y=act_better_FASHOS,
                                                color=texture)) +
  geom_point() +
  geom_errorbar(aes(ymin=act_better_FASHOS-2*act_better_FASHOS_sd,
                    ymax=act_better_FASHOS+2*act_better_FASHOS_sd)) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high))


