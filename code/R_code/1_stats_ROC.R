library(tidyr)
library(dplyr)
library(glmnet)
library(matrixStats)
library(ggplot2)
library(MASS)
source("./analysis_functions.R")

set.seed(2691)

#dataFile <- "../../data/texture_stats/texture_stats.csv"
#dataFile <- "../../data/texture_stats/texture_stats_pixNorm.csv"
dataFile <- "../../data/BSD_stats/BSD_stats_Corr.csv"
saveResults <- "../../data/BSD_results/BSD_angles.Rds"
plotDir <- "../../data/plots/"

# load data
segmentStats <- read.csv(dataFile, sep = ",") %>%
  as_tibble(.)

# get the indices of the different types of stats to use
parNames <- names(segmentStats)
statisticsNames <- get_statistics_names(parNames)
namesPix <- statisticsNames$pixel
namesFAS <- statisticsNames$FAS
namesHOS <- statisticsNames$HOS
designNames <- statisticsNames$design

# normalize the stats
normalizedStats <- dplyr::select(segmentStats, -all_of(designNames)) %>%
  normalize_data(.) %>%
  as_tibble(.)

segmentStatsNorm <- dplyr::select(segmentStats, all_of(designNames)) %>%
  cbind(., normalizedStats) %>%
  as_tibble(.)

# get the angles of the FAS statistics
statsFAS <- dplyr::select(segmentStatsNorm, all_of(c(designNames, namesFAS))) %>%
  compute_angles(.) %>%
  dplyr::mutate(., stats = "FAS", pair = 1:n()) %>%
  as_tibble(.)

statsHOS <- dplyr::select(segmentStatsNorm, all_of(c(designNames, namesHOS))) %>%
  compute_angles(.) %>%
  dplyr::mutate(., stats = "HOS", pair = 1:n()) %>%
  as_tibble(.)

allStats <- rbind(statsFAS, statsHOS) %>%
  pivot_wider(., names_from = stats, values_from = angle)

# calculate the correlations
allCorr <- cor.test(allStats$FAS, allStats$HOS, method = "pearson")
diffCorr <- dplyr::filter(allStats, Type == "diff") %>%
  with(., cor.test(FAS, HOS, method = "pearson"))
sameCorr <- dplyr::filter(allStats, Type == "same") %>%
  with(., cor.test(FAS, HOS, method = "pearson"))

anglePlot <- ggplot(allStats, aes(x = FAS, y = HOS, color = Type)) +
#  geom_abline() +
  geom_point(size = 0.8, alpha = 0.5) +
  theme_bw() +
  labs(x = "Spectral stats angle (rad)", y = "HOS angle (rad)", color = "Segments:") +
  scale_color_manual(labels = c("Different", "Same"),
                       values = c("#b33018", "#14b74b")) +
  xlim(0, 2.6) +
  ylim(0, 2.6) +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid"),
        legend.position = "top", legend.title = element_text(size = 10)) +
  NULL

# discriminant analysis

statsLDA <- lda(Type ~ FAS + HOS, data = allStats)
ldaCoefs <- coef(statsLDA)

randomHOS <- allStats
randomHOS$HOS <- sample(randomHOS$HOS)
randomLDA <- lda(Type ~ FAS + HOS, data = randomHOS)
randomCoefs <- coef(randomLDA)

fischerStats <- dplyr::mutate(allStats, LDA = FAS * ldaCoefs[1] + HOS * ldaCoefs[2])
fischerStats$control <- randomHOS$FAS * randomCoefs[1] + randomHOS$HOS * randomCoefs[2]

rocLDA <- pROC::roc(Type ~ LDA, data = fischerStats)
rocFAS <- pROC::roc(Type ~ FAS, data = fischerStats)
rocHOS <- pROC::roc(Type ~ HOS, data = fischerStats)
rocCtrl <- pROC::roc(Type ~ control, data = fischerStats)

AUC <- list(FAS = rocFAS$auc, HOS = rocHOS$auc, LDA = rocLDA$auc,
            Ctrl = rocCtrl$auc)

coordsFAS <- pROC::coords(rocFAS) %>%
  dplyr::mutate(., type = "FAS")
coordsLDA <- pROC::coords(rocLDA) %>%
  dplyr::mutate(., type = "LDA")
coordsHOS <- pROC::coords(rocHOS) %>%
  dplyr::mutate(., type = "HOS")
coordsCtrl <- pROC::coords(rocCtrl) %>%
  dplyr::mutate(., type = "Ctrl")
coordsAll <- rbind(coordsLDA, coordsFAS, coordsHOS, coordsCtrl)

rocPlot <- dplyr::filter(coordsAll, type != "Ctrl") %>%
  droplevels(.) %>%
  ggplot(., aes(x = specificity, y = sensitivity, color = type)) +
  geom_line() +
  scale_x_reverse(expand = c(0,0)) +
  theme_bw() +
  geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), color = "black") +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_color_manual(name = element_blank(), labels = c("Spectral", "HOS", "Spectral + HOS"),
                       values = c("#2691d4", "#b400e5", "#ba2229")) +
  xlab("Specificity") +
  ylab("Sensitivity") +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid"),
        axis.line.y = element_line(size=0.5, linetype="solid"),
        legend.position = "top", legend.text = element_text(size=8)) +
  NULL

# save results
saveRDS(allStats, saveResults)
ggsave(paste(plotDir, "BSD_angles.png", sep = ""), anglePlot, width = 8,
       height = 9, units = "cm")
ggsave(paste(plotDir, "BSD_roc.png", sep = ""), rocPlot, width = 8.5,
       height = 9, units = "cm")


