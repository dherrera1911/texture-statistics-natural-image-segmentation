library(tidyr)
library(dplyr)
library(glmnet)
library(caret)
library(vcd)
source("./analysis_functions.R")
set.seed(2691)

# data files
dataFileHOS <- "../../data/BSD_stats/analysis_better_HOS_PCA.RDS"
dataFileFASHOS <- "../../data/BSD_stats/analysis_better_FASHOS_PCA.RDS"
dataFileDiff <- "../../data/BSD_stats/analysis_better_FASHOS_Diff_PCA.RDS"

# load data files

############################################################
########## Plot coefficients of PS statistics ##############
############################################################

paramTypeNames <- c("mm", "acr", "acm", "cmc", "pmc", "prc")

coefBetterHOS <- statisticsHOSBetter$coefBetterHOS
coefDf <- data.frame(paramName = rownames(coefBetterHOS[[1]]))
for (l in c(1:length(coefBetterHOS))) {
  colName <- paste("rep", l, sep="")
  coefDf[[colName]] <- as.numeric(coefBetterHOS[[l]])
}
repNamesCol <- paste("rep", c(1:length(coefBetterHOS)), sep="")

coefDf$colMeans <- dplyr::select(coefDf, repNamesCol) %>%
  as.matrix(.) %>%
  rowMeans(.)

coefDf <- dplyr::as_tibble(coefDf)

nonParam <- which(!(coefDf$paramName %in% c(namesFAS, namesHOS)))
coefDf <- coefDf[-nonParam,]

coefDf$paramType <- subset_statistics(paramTypeNames = paramTypeNames,
                               paramVec = coefDf$paramName)

plotDf <- tidyr::pivot_longer(coefDf, cols = repNamesCol,
                              names_to = "rep",
                              values_to = "coefVal")

coefPlot <- ggplot(plotDf, aes(x = factor(paramType), y = abs(coefVal),
                               fill = factor(paramType))) +
  geom_violin() + 
  stat_summary(fun="median", geom="point")


statsSummary <- group_by(plotDf, paramName, paramType) %>%
  summarize(., coefMean = mean(coefVal), coefSD = sd(coefVal)) %>%
  ungroup(.)

coefPlot <- ggplot(statsSummary, aes(x = factor(paramType), y = coefMean,
                               fill = factor(paramType))) +
  geom_violin() +
  stat_summary(fun="median", geom="point")

