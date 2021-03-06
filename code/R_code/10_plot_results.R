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

textureResultsFile <- "../../data/texture_results/2_texture_results.RDS"
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

texturePlotName <- paste(plottingDir, "2_texturePlot.png", sep="")
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

texturePlotName <- paste(plottingDir, "2_texturePlot_noPix.png", sep="")
ggsave(texturePlotName, texturePlot2, width = 10, height = 7, units = "cm") 


##########################
##########################
##### 3 Texture DNN ######
##########################
##########################
textureResultsDNNFile <- "../../data/texture_results/5_texture_dnn.Rds"
textureDNNRes <- readRDS(textureResultsDNNFile) %>%
  dplyr::mutate(., statsName=assign_stat_name_plot(pixel, FAS, HOS),
                error=1-performance)
textureDNNRes$statsName <- factor(textureDNNRes$statsName,
                               levels=c("Pix", "Pix-Spectral",
                                        "Pix-HOS", "Pix-Spectral-HOS",
                                        "Spectral", "HOS", "Spectral-HOS"))

textureDNNPerformanceSummary <- group_by(textureDNNRes, statsName, pixel) %>%
  summarize(., meanError=mean(error), sdError=sd(error))

#### plot results when using pixel data ####
textureDNNPlot <- dplyr::filter(textureDNNRes, pixel==1) %>%
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
  scale_x_discrete(name="Parameters") +
  NULL

textureDNNPlotName <- paste(plottingDir, "5_textureDNNPlot.png", sep="")
ggsave(textureDNNPlotName, textureDNNPlot, width=10, height=7, units = "cm") 

#### plot results without pixel data ####
textureDNNPlot2 <- dplyr::filter(textureDNNRes, pixel==0) %>%
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
  scale_x_discrete(name="Parameters") +
  NULL

texturePlotName <- paste(plottingDir, "5_texturePlot_noPix.png", sep="")
ggsave(texturePlotName, textureDNNPlot2, width=10, height=7, units = "cm") 


#############################
#############################
##### 4 BSD results #########
#############################
#############################

bsdFile <- "../../data/BSD_results/3_BSD_results_noNorm.RDS"
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
  ylab("Error rate (%)") +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid")) +
  scale_x_discrete(name = "Parameters")

BSDPlotName <- paste(plottingDir, "BSDPlot.png", sep="")
ggsave(BSDPlotName, bsdPlot, width=10, height=7, units = "cm") 

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
  ylab("Error rate (%)") +
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
  summarize(., meanError=mean(error), sdError=sd(error), count=n())

bsdDNNPlot <- dplyr::filter(bsdDNNRes, pixel==1) %>%
  droplevels(.) %>%
  ggplot(aes(x=statsName, y=error*100, color=statsName)) +
  geom_jitter(width=0.2, height=0, alpha=0.5, shape=1) +
  scale_color_manual(name="Statistic groups",
                       values = c("#000000", "#2691d4", "#b400e5", "#ba2229"),
                     guide = "none") +
  stat_summary(fun = "mean") +
  ylim(0, 30) +
  theme_bw() +
  theme(panel.border=element_blank(),
        axis.line.x=element_line(size=0.5, linetype="solid"),
        axis.line.y=element_line(size=0.5, linetype="solid")) +
  ylab("Error rate (%)") +
  scale_x_discrete(name="Parameters")

bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot.png", sep="")
ggsave(bsdDNNPlotName, bsdDNNPlot, width=10, height=7, units="cm") 

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
  ylab("Error rate (%)") +
  scale_x_discrete(name = "Parameters")

bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot_noPix.png", sep="")
ggsave(bsdDNNPlotName, bsdDNNPlot2, width=10, height=7, units="cm") 

##################################
##################################
##### 5.2 BSD DNN params #########
##################################
##################################

#bsdDNNFile <- "../../data/BSD_results/4_BSD_dnn_params.Rds"
#bsdDNNRes <- readRDS(bsdDNNFile) %>%
#  dplyr::mutate(., statsName=assign_stat_name_plot(pixel, FAS, HOS),
#                error=1-performance)
#bsdDNNRes$statsName <- factor(bsdDNNRes$statsName,
#                               levels=c("Pix", "Pix-Spectral",
#                                        "Pix-HOS", "Pix-Spectral-HOS",
#                                        "Spectral", "HOS", "Spectral-HOS"))
#
#bsdDNNPerformanceSummary <- group_by(bsdDNNRes, statsName, pixel,
#                                     architecture, regularization, epochs) %>%
#  summarize(., meanError=mean(error), sdError=sd(error), count=n())
#
## compare performance of NN to linear model
#linearSpectralHOS <- dplyr::filter(bsdPerformanceSummary,
#                                   statsName=="Pix-Spectral-HOS")[["meanError"]]
#
#performanceDiff <- dplyr::filter(bsdDNNPerformanceSummary,
#                                 statsName=="Pix-Spectral-HOS") %>%
#  dplyr::mutate(., errorDiff=linearSpectralHOS-meanError)
#
#bestFASHOS <- performanceDiff[which.min(performanceDiff$errorDiff),]
#
#
##### compare NN regimes
#allEpochs <- unique(bsdDNNRes$epochs)
#allReguls <- unique(bsdDNNRes$regularization)
#bsdDNNParPlot <- list()
#bsdDNNParPlot2 <- list()
#for (ep in c(1:length(allEpochs))) {
#  bsdDNNParPlot[[ep]] <- list()
#  bsdDNNParPlot2[[ep]] <- list()
#  for (reg in c(1:length(allReguls))) {
#    bsdDNNParPlot[[ep]][[reg]] <- dplyr::filter(bsdDNNRes, pixel==1,
#                                       epochs==allEpochs[ep] &
#                                         regularization==allReguls[reg]) %>%
#      droplevels(.) %>%
#      ggplot(aes(x=statsName, y=error*100, color=statsName)) +
#      geom_jitter(width=0.2, height=0, alpha=0.5, shape=1) +
#      scale_color_manual(name = "Statistic groups",
#                           values = c("#000000", "#2691d4", "#b400e5", "#ba2229"),
#                         guide = "none") +
#      stat_summary(fun = "mean") +
#      ylim(0, 30) +
#      theme_bw() +
#      facet_wrap(facets=c("architecture", "regularization")) +
#      theme(panel.border = element_blank(),
#            axis.line.x = element_line(size=0.5, linetype="solid"),
#            axis.line.y = element_line(size=0.5, linetype="solid")) +
#      ylab("Error rate (%)") +
#      scale_x_discrete(name = "Parameters",
#                       labels = c("Pix", "Pix-Spectral", "Pix-HOS",
#                                  "Pix-Spectral-HOS"))
#    bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot_param_", 
#                            allEpochs[ep], "_", allReguls[reg],
#                            ".png", sep="")
#    ggsave(bsdDNNPlotName, bsdDNNParPlot[[ep]][[reg]],
#           width=20, height=14, units = "cm") 
#    bsdDNNParPlot2[[ep]][[reg]] <- dplyr::filter(bsdDNNRes, pixel==0,
#                                       epochs==allEpochs[ep] &
#                                         regularization==allReguls[reg]) %>%
#      droplevels(.) %>%
#      ggplot(aes(x=statsName, y=error*100, color=statsName)) +
#      geom_jitter(width=0.2, height=0, alpha=0.5, shape=1) +
#      scale_color_manual(name = "Statistic groups",
#                           values = c("#2691d4", "#b400e5", "#ba2229"),
#                         guide = "none") +
#      stat_summary(fun = "mean") +
#      ylim(0, 30) +
#      theme_bw() +
#      facet_wrap(facets=c("architecture", "regularization")) +
#      theme(panel.border = element_blank(),
#            axis.line.x = element_line(size=0.5, linetype="solid"),
#            axis.line.y = element_line(size=0.5, linetype="solid")) +
#      ylab("Error rate (%)") +
#      scale_x_discrete(name = "Parameters")
#    bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot_param_", 
#                            allEpochs[ep], "_", allReguls[reg],
#                            "_noPix.png", sep="")
#    ggsave(bsdDNNPlotName, bsdDNNParPlot2[[ep]][[reg]],
#           width=20, height=14, units="cm") 
#  }
#}

###################################################
###################################################
##### 5.3 BSD DNN convergence progression #########
###################################################
###################################################
dnnHistFile <- "../../data/BSD_results/4_BSD_dnn_params_history.Rds"
dnnHistoryList <- readRDS(dnnHistFile)

history1 <- readRDS("../../data/BSD_results/4_BSD_dnn_params_history_1.Rds")
history2 <- readRDS("../../data/BSD_results/4_BSD_dnn_params_history_2.Rds")
history3 <- readRDS("../../data/BSD_results/4_BSD_dnn_params_history_3_reg01.Rds")
history4 <- readRDS("../../data/BSD_results/4_BSD_dnn_params_history_4_reg01.Rds")
history5 <- readRDS("../../data/BSD_results/4_BSD_dnn_params_history_5_reg001.Rds")
history6 <- readRDS("../../data/BSD_results/4_BSD_dnn_params_history_6_reg01_stats.Rds")

historyList <- list(history1, history2, history3, history4, history5, history6)
historyDfList <- list()
for (l in c(1:length(historyList))) {
  architectures <- names(historyList[[l]])
  statsNames <- names(historyList[[l]][[1]])
  historyDf <- data.frame()
  for (arch in architectures) {
    for (sn in statsNames) {
      errorMeans <- 1-(colMeans(historyList[[l]][[arch]][[sn]]))
      errorSd <- sqrt(colVars(historyList[[l]][[arch]][[sn]]))
      tempDf <- data.frame(architecture=arch, stats=sn,
                           epoch=c(1:length(errorMeans)),
                           error=errorMeans, sd=errorSd)
      historyDf <- rbind(historyDf, tempDf)
    }
  }
  historyDfList[[l]] <- historyDf
}
# add the regularization of each run
historyDfList[[1]]$reg <- 0.003
historyDfList[[2]]$reg <- 0.003
historyDfList[[3]]$reg <- 0.01
historyDfList[[4]]$reg <- 0.01
historyDfList[[5]]$reg <- 0.001
# bind all together
historyDf <- rbind(historyDfList[[1]], historyDfList[[2]],
                   historyDfList[[3]], historyDfList[[4]],
                   historyDfList[[5]])

# Plot the performance progression
linearHOS <- dplyr::filter(bsdPerformanceSummary,
                           statsName=="Pix-HOS")[["meanError"]]
historyPlot_HOS <- dplyr::filter(historyDf, stats=="Pix_HOS" & reg==0.001) %>%
  ggplot(., aes(x=epoch, y=error*100, color=architecture)) +
  #facet_grid(~reg) +
  geom_line() +
  ylab("Validation error rate (%)") +
  xlab("Epoch") +
  theme_bw() +
  ylim(c(0, 60)) +
  labs(color="Architecture") +
  geom_hline(yintercept=linearHOS*100, linetype=2)
dnnHistoryPlotName <- paste(plottingDir, "4_NNhistoryPlot_HOS.png", sep="")
ggsave(dnnHistoryPlotName, historyPlot_HOS, width=12, height=9, units="cm") 

# Plot the performance progression
linearFASHOS <- dplyr::filter(bsdPerformanceSummary,
                           statsName=="Pix-Spectral-HOS")[["meanError"]]
historyPlot_FASHOS <- dplyr::filter(historyDf, stats=="Pix_FAS_HOS") %>%
  ggplot(., aes(x=epoch, y=error*100, color=architecture)) +
  facet_grid(~reg) +
  geom_line() +
  ylab("Error") +
  xlab("Epoch") +
  theme_bw() +
  ylim(c(0, 50)) +
  geom_hline(yintercept=linearFASHOS*100, linetype=2)
dnnHistoryPlotName <- paste(plottingDir, "4_NNhistoryPlot_FASHOS.png", sep="")
ggsave(dnnHistoryPlotName, historyPlot_FASHOS, width=10, height=7, units="cm") 

dnnPlotBSD <- ggpubr::ggarrange(plotlist=list(historyPlot_HOS, bsdDNNPlot),
                                nrow=1, common.legend=FALSE,
                                labels=list("(A)", "(B)"))

dnnFigPlotName <- paste(plottingDir, "4_dnn_duplex.png", sep="")
ggsave(dnnFigPlotName, dnnPlotBSD, width=19, height=7.5, units="cm") 

# Plot for one architecture, all the statistics progression
#linearFASHOS <- dplyr::filter(bsdPerformanceSummary, )[["meanError"]]
historyPlot_allStats <- dplyr::filter(historyDfList[[6]]) %>%
  ggplot(., aes(x=epoch, y=error*100, color=stats, linetype=architecture)) +
  #facet_grid(~reg) +
  geom_line() +
  ylab("Error") +
  xlab("Epoch") +
  theme_bw() +
  ylim(c(0, 50)) +
#  geom_hline(yintercept=linearFASHOS*100, linetype=2)
  NULL

dnnHistoryPlotName <- paste(plottingDir, "4_NNhistoryPlot_FASHOS.png", sep="")
ggsave(dnnHistoryPlotName, historyPlot_FASHOS, width=10, height=7, units="cm") 

##################################
##################################
##### 6 analyze agreement ########
##################################
##################################

agreementFile <- "../../data/BSD_results/6_classification_agreement.Rds"
#agreementFile <- "../../data/texture_results/6_classification_agreement_textures.Rds"
agreementDf <- readRDS(agreementFile) %>%
  dplyr::mutate(., errorAgreementHOS=wrongHOS*wrongFAS,
                   errorAgreementFASHOS=wrongFASHOS*wrongHOS)

agreementSel <- dplyr::select(agreementDf, same, betterHOS, betterFASHOS,
                              wrongFAS, errorAgreementHOS, errorAgreementFASHOS,
                              wrongHOS, wrongFASHOS) %>%
  tidyr::pivot_longer(., cols=c("betterHOS", "betterFASHOS", "wrongFAS",
                                "errorAgreementHOS", "errorAgreementFASHOS",
                                "wrongHOS", "wrongFASHOS"),
                      names_to="comparison", values_to="result")

agreementSummaryCoarse1 <- agreementSel %>%
  group_by(., comparison) %>% 
  summarize(., proportion=mean(result))

agreementSummaryCoarse <- agreementSel %>%
  group_by(., comparison) %>% 
  summarize(., proportion=mean(result)) %>%
  dplyr::filter(., comparison %in% c("betterHOS", "betterFASHOS", "wrongFAS") %>%
  dplyr::mutate(., comparison=factor(comparison,
                                       levels=c("betterHOS", "betterFASHOS",
                                                "wrongFAS")
                                       labels=c("HOS", "Spectral-HOS",
                                                "Bad-Spectral")))

agreementSummary <- agreementSel %>%
  group_by(., same, comparison) %>% 
  summarize(., proportion=mean(result)) %>%
  dplyr::mutate(., comparison=factor(comparison,
                                       levels=c("betterHOS", "betterFASHOS",
                                                "wrongFAS"),
                                       labels=c("HOS", "Spectral-HOS",
                                                "Bad-Spectral")),
                same=factor(same, levels=c(0, 1),
                            labels=c("Different", "Same")))

agreementPlot <- dplyr::filter(agreementSummary, outcome=="Better") %>%
  ggplot(data=., aes(x=same, y=proportion*100, fill=same)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~comparison) +
  scale_fill_manual(name="Pair class",
                     labels=c("Different", "Same"),
                       values=c("#b33018", "#14b74b")) +
  scale_x_discrete(labels=NULL) +
  scale_y_continuous(limits=c(0,10), expand=c(0, 0)) +
  xlab("") +
  ylab("% pairs better than spectral") +
  theme_bw() +
  theme(axis.ticks=element_blank())

agreementPlotName <- paste(plottingDir, "agreement_bars.png", sep="")
ggsave(agreementPlotName,  agreementPlot, width=13, height=9, units="cm") 

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

agreementPredTextureFile <- "../../data/texture_results/7_confusion_matrices_texture.RDS"
agreementPredTextureDf <- readRDS(agreementPredTextureFile)

predictionHOSTexture <- caret::confusionMatrix(agreementPredTextureDf$confusionHOS)
predictionFASHOSTexture <- caret::confusionMatrix(agreementPredTextureDf$confusionFASHOS)


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

#subsetsFile <- "../../data/BSD_results/8_subset_stats_performance.RDS"
subsetsFile <- "../../data/texture_results/8_subset_stats_performance_texture.RDS"
subsetsDf <- readRDS(subsetsFile) %>%
  dplyr::mutate(., statsName=assign_subset_name_plot(FAS, acm, cmc, pmc, prc),
                error=1-performance, numberSubsets=FAS+acm+cmc+pmc+prc)

subsetsSummary <- group_by(subsetsDf, statsName, FAS, acm, cmc, pmc, prc, numberSubsets) %>%
  summarize(., meanError=mean(error), sdError=sd(error)/sqrt(n())) %>%
  dplyr::mutate(., numberSubsets=factor(numberSubsets))
subsetsSummary <- subsetsSummary[order(subsetsSummary$numberSubsets),]

spectralPerformance <- dplyr::filter(subsetsSummary, numberSubsets==1 & FAS==1)[["meanError"]]
fullModelPerformance <- dplyr::filter(subsetsSummary, numberSubsets==5)[["meanError"]]

# filter data of interest
singleStats <- dplyr::filter(subsetsSummary, numberSubsets==1 & FAS==0)
addedStats <- dplyr::filter(subsetsSummary, numberSubsets==2 & FAS==1)
removedStats <- dplyr::filter(subsetsSummary, numberSubsets==4 & FAS==1)
# put friendlier names
addedStats$statsName <- sub("Spectral_", "+ ", addedStats$statsName)
statsNames <- c("Space", "Orientation", "Scale", "Phase")
indVec <- NULL
nameVec <- NA
for (sN in statsNames) {
  ind <- which(!(c(1:4) %in% grep(sN, removedStats$statsName)))
  nameVec[ind] <- sN
}
removedStats$statsName <- paste("-", nameVec)

#plotLims <- c(0, 38)
#plotLims <- c(0, 31)
singleStatsPlot <- singleStats %>%
  ggplot(data=., aes(x=statsName, y=meanError*100)) +
  geom_bar(stat="identity", width=0.5) +
#  geom_errorbar(aes(ymin=(meanError-sdError*2)*100,
#                    ymax=(meanError+sdError*2)*100), width=0.1) +
  geom_hline(yintercept=spectralPerformance*100, color="#2691d4", linetype=2) +
  geom_hline(yintercept=fullModelPerformance*100, color="#ba2229", linetype=2) +
  ylab("Error rate (%)") +
  xlab("Statistics used") +
#  ylim(plotLims) +
  ggtitle("Subsets alone") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle=30, vjust=0.8, hjust=0.8),
        axis.title.x = element_blank(),
        title = element_text(size=8)) +
  xlab("") +
  NULL

addedStatsPlot <- addedStats %>%
  ggplot(data=., aes(x=statsName, y=meanError*100)) +
  geom_bar(stat="identity", width = 0.5) +
#  geom_errorbar(aes(ymin=(meanError-sdError*2)*100,
#                    ymax=(meanError+sdError*2)*100), width=0.1) +
  geom_hline(yintercept=spectralPerformance*100, color="#2691d4", linetype=2) +
  geom_hline(yintercept=fullModelPerformance*100, color="#ba2229", linetype=2) +
  ylab("Error rate (%)") +
  xlab("Statistics used") +
  #ylim(plotLims) +
  ggtitle("Spectral +subsets") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=30, vjust=0.8, hjust=0.8),
        axis.title.x = element_blank(),
        title = element_text(size=8)) +
  xlab("") +
  NULL

removedStatsPlot <- removedStats %>%
  ggplot(data=., aes(x=statsName, y=meanError*100)) +
  geom_bar(stat="identity", width = 0.5) +
#  geom_errorbar(aes(ymin=(meanError-sdError*2)*100,
#                    ymax=(meanError+sdError*2)*100), width=0.1) +
  geom_hline(yintercept=spectralPerformance*100, color="#2691d4", linetype=2) +
  geom_hline(yintercept=fullModelPerformance*100, color="#ba2229", linetype=2) +
  ylab("Error rate (%)") +
  xlab("Statistics used") +
  #ylim(plotLims) +
  ggtitle("Full model -subsets") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle=30, vjust=0.8, hjust=0.8),
        axis.title.x = element_blank(),
        title = element_text(size=8)) +
  xlab("") +
  NULL

allPlots <- ggpubr::ggarrange(plotlist=list(singleStatsPlot, addedStatsPlot,
                                            removedStatsPlot), nrow=1, common.legend=FALSE,
  labels=list("(A)", "(B)", "(C)"))

subsetsPlotName <- paste(plottingDir, "8_subsets_BSD.png", sep="")
#subsetsPlotName <- paste(plottingDir, "8_subsets_texture.png", sep="")
ggsave(subsetsPlotName, allPlots, width=16, height=6, units="cm") 

# pairwise subsets
singlesSummary <- dplyr::filter(subsetsSummary, numberSubsets==1 & FAS==0)
pairsSummary <- dplyr::filter(subsetsSummary, numberSubsets==2 & FAS==0)

######################################
##### 9 Experiment results  ##########
######################################

expResultsFile <- "../../data/experiment_stats/experimentResults.RDS"
expSegmentationFile <- "../../data/BSD_results/9_segmentation_predictions.RDS"
expAgreementFile <- "../../data/BSD_results/9_better_stats_predictions.RDS"
#expSegmentationFile <- "../../data/texture_results/9_segmentation_predictions_texture.RDS"
#expAgreementFile <- "../../data/texture_results/9_better_stats_predictions_texture.RDS"

expResults <- readRDS(expResultsFile) %>%
  as_tibble(.)
expSeg <- readRDS(expSegmentationFile) %>%
  as_tibble(.)

summaryExpSeg <- group_by(expSeg, texture) %>%
  summarize(., mean_same_HOS=mean(same_HOS),
            mean_same_FASHOS=mean(same_FASHOS),
            act_diff_HOS=-mean(activation_HOS),
            act_diff_FASHOS=-mean(activation_FASHOS),
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
  dplyr::mutate(., estimate=-estimate, conf.low=-conf.low, conf.high=-conf.high) %>%
  dplyr::rename(., Texture=texture)
  
                #texture=factor(texture, labels=c("T1", "T2", "T3", "T4")),
                

experimentHOSSegPlot <- ggplot(expDataFull, aes(x=estimate, y=act_diff_HOS,
                                                color=Texture)) +
  geom_point() +
  geom_errorbar(aes(ymin=act_diff_HOS-2*act_diff_HOS_sd,
                    ymax=act_diff_HOS+2*act_diff_HOS_sd), size=1, width=0) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), size=1, height=0) +
  xlab("Experimental HOS segmentation") +
  ylab("Model segmentation (linear activation)") +
  theme_bw() +
  NULL

experimentFASHOSSegPlot <- ggplot(expDataFull, aes(x=estimate, y=act_diff_FASHOS,
                                                color=Texture)) +
  geom_point() +
  geom_errorbar(aes(ymin=act_diff_FASHOS-2*act_diff_FASHOS_sd,
                    ymax=act_diff_FASHOS+2*act_diff_FASHOS_sd), size=1, width=0) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), size=1, height=0) +
  xlab("Experimental HOS segmentation") +
  ylab("Model segmentation (linear activation)") +
  theme_bw() +
  NULL

experimentHOSAgreementPlot <- ggplot(expDataFull, aes(x=estimate, y=act_better_HOS,
                                                color=Texture)) +
  geom_point() +
  geom_errorbar(aes(ymin=act_better_HOS-2*act_better_HOS_sd,
                    ymax=act_better_HOS+2*act_better_HOS_sd), size=1, width=0) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), size=1, height=0) +
  xlab("Experimental HOS segmentation") +
  ylab("HOS usefulness (linear activation)") +
  theme_bw()

experimentFASHOSAgreementPlot <- ggplot(expDataFull, aes(x=estimate, y=act_better_FASHOS,
                                                color=Texture)) +
  geom_point() +
  geom_errorbar(aes(ymin=act_better_FASHOS-2*act_better_FASHOS_sd,
                    ymax=act_better_FASHOS+2*act_better_FASHOS_sd), size=1, width=0) +
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), size=1, height=0) +
  xlab("Experimental HOS segmentation") +
  ylab("HOS usefulness (linear activation)") +
  theme_bw()

experimentHOSPlot <- ggpubr::ggarrange(plotlist = list(experimentHOSSegPlot,
                                            experimentHOSAgreementPlot),
                                       nrow=1, common.legend=TRUE,
                                       labels=list("(A)", "(B)"), legend="right")

experimentFASHOSPlot <- ggpubr::ggarrange(plotlist = list(experimentFASHOSSegPlot,
                                            experimentFASHOSAgreementPlot),
                                       nrow=1, common.legend=TRUE,
                                       labels=list("(A)", "(B)"), legend="right")

saveHOSName <- paste("../../data/plots/9_experiment_model", "_HOS", ".png", sep="")
#saveHOSName <- paste("../../data/plots/9_experiment_model", "_HOS", "_texture.png", sep="")
ggsave(saveHOSName, experimentHOSPlot, width=23, height=10, units="cm") 
saveFASHOSName <- paste("../../data/plots/9_experiment_model", "_FASHOS", ".png", sep="")
#saveFASHOSName <- paste("../../data/plots/9_experiment_model", "_FASHOS", "_texture.png", sep="")
ggsave(saveFASHOSName, experimentFASHOSPlot, width=22, height=10, units="cm") 


