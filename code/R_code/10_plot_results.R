library(tidyr)
library(dplyr)
library(ggplot2)
library(matrixStats)
source("./analysis_functions.R")
set.seed(2691)
#library(Hmisc)

plottingDir <- "../../data/plots/"
modelOrder <- c("Pix", "FA", "HOS", "FA_HOS")

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

parPlot <- ggplot(data = nParDf, aes(x = Statistics, y = nPar, fill = Statistics)) +
  geom_bar(stat="identity", width = 0.5) +
  scale_fill_manual(name = "Statistic groups",
                     labels = c("Pixel", "FAS", "HOS"),
                       values = c("#000000", "#2691d4", "#ba2229"),
                     guide = "none") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid")) +
  ylab("Number of parameters") +
  scale_x_discrete(name = "Parameters",
                   labels = c("Pixel", "FAS", "HOS"))

ggsave(paste(plottingDir, "nParameters.png", sep=""), parPlot,
       width = 10, height = 10, units = "cm") 

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
  scale_x_discrete(name = "Parameters",
                   labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS")) +
  NULL

texturePlotName <- paste(plottingDir, "texturePlot.png", sep="")
ggsave(texturePlotName, texturePlot, width = 10, height = 7, units = "cm") 


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

#### plot results when using pixel data ####
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
  scale_x_discrete(name = "Parameters",
                   labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS"))

BSDPlotName <- paste(plottingDir, "BSDPlot.png", sep="")
ggsave(BSDPlotName, bsdPlot, width = 10, height = 7, units = "cm") 

#########################
#########################
##### 5 BSD DNN #########
#########################
#########################

bsdDNNFile <- "../../data/BSD_results/4_BSD_dnn.Rds"
bsdDNNRes <- readRDS(bsdDNNFile) %>%
  dplyr::mutate(., statsName=assign_stat_name_plot(pixel, FAS, HOS),
                error=1-performance)
bsdDNNRes$statsName <- factor(bsdDNNRes$statsName,
                               levels=c("Pix", "Pix-Spectral",
                                        "Pix-HOS", "Pix-Spectral-HOS",
                                        "Spectral", "HOS", "Spectral-HOS"))

bsdDNNPerformanceSummary <- group_by(bsdDNNRes, statsName, pixel,
                                     architecture, regularization, epochs) %>%
  summarize(., meanError=mean(error), sdError=sd(error),
            minrep=max(rep))

# compare performance of NN to linear model
linearSpectralHOS <- dplyr::filter(bsdPerformanceSummary,
                                   statsName=="Pix-Spectral-HOS")[["meanError"]]

performanceDiff <- dplyr::filter(bsdDNNPerformanceSummary,
                                 statsName=="Pix-Spectral-HOS") %>%
  dplyr::mutate(., errorDiff=linearSpectralHOS-meanError)

#### compare NN regimes
allEpochs <- unique(bsdDNNRes$epochs)
bsdDNNPlots <- list()
for (ep in c(1:length(allEpochs))) {
  bsdDNNPlots[[ep]] <- dplyr::filter(bsdDNNRes, pixel==1,
                                     epochs==allEpochs[ep]) %>%
    droplevels(.) %>%
    ggplot(aes(x = statsName, y = error  * 100, color = statsName)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5, shape = 1) +
    scale_color_manual(name = "Statistic groups",
                       labels = c("Pix", "Pix-Spectral", "Pix-HOS",
                                  "Pix-Spectral-HOS"),
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
}


#### plot main results ####
bsdDNNPlot <- dplyr::filter(bsdDNNRes, pixel==1) %>%
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
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid")) +
  ylab("error rate (%)") +
  scale_x_discrete(name = "Parameters",
                   labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS"))


bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot.png", sep="")
ggsave(bsdDNNPlotName, bsdDNNPlot, width = 10, height = 7, units = "cm") 


##################################
##### 6 analyze agreement #########
##################################
agreementFile <- "../../data/BSD_results/6_classification_agreement.csv"
agreementDf <- read.csv(agreementFile, stringsAsFactors=FALSE)

agreementSummary <- summarize(agreementDf, betterHOS=mean(betterHOS),
                               betterFASHOS=mean(betterFASHOS))


##################################
##### 7 predict agreement ########
##################################
agreementFile <- "../../data/BSD_results/8_subset_stats_performance.RDS"
agreementDf <- readRDS(agreementFile)

agreementSummary <- summarize(agreementDf, betterHOS=mean(betterHOS),
                               betterFASHOS=mean(betterFASHOS))




