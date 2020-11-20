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
##### Plot N of parameters ########
##################################
statsFile <- "../../data/texture_stats/texture_stats_statsNorm.csv"
statsDf <- read.csv(statsFile, sep = ",") %>%
  as_tibble(.) %>%
  dplyr::select(., -texture, -quadrant)

parNames <- names(statsDf)
noVarInds <- which(colVars(as.matrix(statsDf)) == 0)
noVarNames <- parNames[noVarInds]
parNames <- parNames[-noVarInds]
indPix <- c(grep("pix", parNames), grep("^LP", parNames))
indFAS <- c(grep("mm", parNames), grep("acr", parNames))
indHOS <- c(grep("acm", parNames), grep("cmc", parNames),
             grep("pmc", parNames), grep("prc", parNames))
nPix <- length(indPix)
nFAS <- length(indFAS)
nHOS <- length(indHOS)

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
##### Texture results ########
#############################

textureResultsFile <- "../../data/texture_results/texture_results.Rds"

textureRes <- readRDS(textureResultsFile) %>%
  pivot_longer(., c("Pix", "FA", "HOS", "FA_HOS"), names_to = "Model_params",
               values_to = "Performance") %>%
  dplyr::filter(., !is.na(Performance))

# add clearer the stats
colFAS <- rep(0, nrow(textureRes))
colFAS[grep("FA", textureRes$Model_params)] <- 1
colHOS <- rep(0, nrow(textureRes))
colHOS[grep("*HOS", textureRes$Model_params)] <- 1

textureRes$FAS <- colFAS
textureRes$HOS <- colHOS

# tidying
pixRes <- dplyr::filter(textureRes, usePix == 1) %>%
  dplyr::mutate(., Model_params =
                c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS")[usePix + FAS + 2*HOS])

modelOrder <- c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS", "FA", "HOS", "FA_HOS")

textureRes <- dplyr::filter(textureRes, usePix == 0) %>%
  rbind(., pixRes) %>%
  dplyr::mutate(., Error = 1 - Performance)

textureRes$Model_params <- factor(textureRes$Model_params, levels = modelOrder)

# summarize performance
texturePerformance <- textureRes %>%
  dplyr::select(., Model_params, Error) %>%
  group_by(., Model_params) %>%
  dplyr::summarize(., mean_error = mean(Error), error_sd = sd(Error))


#### plot results when using pixel data ####
texturePlot <- dplyr::filter(textureRes, usePix == 1) %>%
  ggplot(aes(x = Model_params, y = Error  * 100, color = Model_params)) +
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


#############################
##### Texture DNN ########
#############################

textureResultsDNNFile <- "../../data/texture_results/texture_dnn.Rds"

textureDNNRes <- readRDS(textureResultsDNNFile) %>%
  pivot_longer(., c("Pix", "FA", "HOS", "FA_HOS"), names_to = "Model_params",
               values_to = "Performance") %>%
  dplyr::filter(., !is.na(Performance))

# add clearer the stats
colFAS <- rep(0, nrow(textureDNNRes))
colFAS[grep("FA", textureDNNRes$Model_params)] <- 1
colHOS <- rep(0, nrow(textureDNNRes))
colHOS[grep("*HOS", textureDNNRes$Model_params)] <- 1

textureDNNRes$FAS <- colFAS
textureDNNRes$HOS <- colHOS

pixRes <- dplyr::filter(textureDNNRes, usePix == 1) %>%
  dplyr::mutate(., Model_params =
                c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS")[usePix + FAS + 2*HOS])

modelOrder <- c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS", "FA", "HOS", "FA_HOS")

textureDNNRes <- dplyr::filter(textureDNNRes, usePix == 0) %>%
  rbind(., pixRes) %>%
  dplyr::mutate(., Error = 1 - Performance)

textureDNNRes$Model_params <- factor(textureDNNRes$Model_params, levels = modelOrder)

texturePerformanceDNN <- textureDNNRes %>%
  dplyr::select(., Model_params, Error) %>%
  group_by(., Model_params) %>%
  dplyr::summarize(., mean_error = mean(Error), error_sd = sd(Error))


#### plot results when using pixel data ####
textureDNNPlot <- dplyr::filter(textureDNNRes, usePix == 1) %>%
  ggplot(aes(x = Model_params, y = Error  * 100, color = Model_params)) +
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
  ylab("Error rate (%)") +
  scale_x_discrete(name = "Parameters",
                   labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS")) +
  NULL

textureDNNPlotName <- paste(plottingDir, "textureDNNPlot.png", sep="")
ggsave(textureDNNPlotName, textureDNNPlot, width = 10, height = 7, units = "cm") 


#############################
##### BSD results ########
#############################

bsdFile <- "../../data/BSD_results/BSD_results.Rds"
modelOrder <- c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS", "FA", "HOS", "FA_HOS")
bsdRes <- readRDS(bsdFile) %>%
  pivot_longer(., c("Pix", "FA", "HOS", "FA_HOS"), names_to = "Model_params",
               values_to = "Performance") %>%
  dplyr::filter(., !is.na(Performance))

colFAS <- rep(0, nrow(bsdRes))
colFAS[grep("FA", bsdRes$Model_params)] <- 1
colHOS <- rep(0, nrow(bsdRes))
colHOS[grep("*HOS", bsdRes$Model_params)] <- 1

bsdRes$FAS <- colFAS
bsdRes$HOS <- colHOS

pixRes <- dplyr::filter(bsdRes, usePix == 1) %>%
  dplyr::mutate(., Model_params =
                c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS")[usePix + FAS + HOS*2])

bsdRes <- dplyr::filter(bsdRes, usePix == 0) %>%
  rbind(., pixRes) %>%
  dplyr::mutate(., Error = 1 - Performance)

# order the factors
bsdRes$Model_params <- factor(bsdRes$Model_params, levels = modelOrder)

bsdPerformance <- bsdRes %>%
  dplyr::select(., Model_params, Error) %>%
  group_by(., Model_params) %>%
  dplyr::summarize(., mean_error = mean(Error), error_sd = sd(Error))

#### plot results when using pixel data ####
bsdPlot <- bsdRes %>%
  dplyr::filter(., usePix == 1) %>%
  droplevels(.) %>%
  ggplot(aes(x = Model_params, y = Error  * 100, color = Model_params)) +
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
  scale_x_discrete(name = "Parameters",
                   labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS"))

BSDPlotName <- paste(plottingDir, "BSDPlot.png", sep="")
ggsave(BSDPlotName, bsdPlot, width = 10, height = 7, units = "cm") 


#############################
##### BSD DNN ########
#############################

bsdDNNFile <- "../../data/BSD_results/BSD_dnn.Rds"
modelOrder <- c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS", "FA", "HOS", "FA_HOS")
bsdDNNRes <- readRDS(bsdDNNFile) %>%
  pivot_longer(., c("Pix", "FA", "HOS", "FA_HOS"), names_to = "Model_params",
               values_to = "Performance") %>%
  dplyr::filter(., !(usePix == 0 & Model_params == "Pix")) %>%
  dplyr::mutate(., Error = 1 - Performance)

colFAS <- rep(0, nrow(bsdDNNRes))
colFAS[grep("FA", bsdDNNRes$Model_params)] <- 1
colHOS <- rep(0, nrow(bsdDNNRes))
colHOS[grep("*HOS", bsdDNNRes$Model_params)] <- 1

bsdDNNRes$FAS <- colFAS
bsdDNNRes$HOS <- colHOS

pixRes <- dplyr::filter(bsdDNNRes, usePix == 1) %>%
  dplyr::mutate(., Model_params =
                c("Pix", "Pix_FAS", "Pix_HOS", "Pix_FAS_HOS")[usePix + FAS + HOS*2])

bsdDNNRes <- dplyr::filter(bsdDNNRes, usePix == 0) %>%
  rbind(., pixRes)

# order the factors
bsdDNNRes$Model_params <- factor(bsdDNNRes$Model_params, levels = modelOrder)

bsdDNNPerformance <- bsdDNNRes %>%
 # dplyr::filter(., usePix == 1) %>%
#  droplevels(.) %>%
  dplyr::select(., Model_params, Error) %>%
  group_by(., Model_params) %>%
  dplyr::summarize(., mean_error = mean(Error), error_sd = sd(Error))


#### plot results ####
bsdDNNPlot <- dplyr::filter(bsdDNNRes, usePix == 1) %>%
  droplevels(.) %>%
  ggplot(aes(x = Model_params, y = Error  * 100, color = Model_params)) +
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
  ylab("Error rate (%)") +
  scale_x_discrete(name = "Parameters",
                   labels = c("Pix", "Pix-Spectral", "Pix-HOS", "Pix-Spectral-HOS"))


bsdDNNPlotName <- paste(plottingDir, "bsdDNNPlot.png", sep="")
ggsave(bsdDNNPlotName, bsdDNNPlot, width = 10, height = 7, units = "cm") 

