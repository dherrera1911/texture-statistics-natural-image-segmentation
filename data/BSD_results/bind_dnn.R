allFiles <- dir()
fileNames <- allFiles[grep("4_BSD", allFiles)]
fileNames <- fileNames[-grep("single", fileNames)]
mainFile <- grep("4_BSD_dnn.R", fileNames)
if (length(mainFile) > 0) {
  fileNames <- fileNames[-mainFile]
}

allDNN <- NULL
for(fn in fileNames) {
  tempDf <- readRDS(fn)
  allDNN <- rbind(allDNN, tempDf)
}
saveRDS(allDNN, "4_BSD_dnn.RDS")

fas8 <- readRDS("./8_subset_stats_performance_FAS.RDS")
nofas8 <- readRDS("./8_subset_stats_performance_noFAS.RDS")

all8 <- rbind(fas8, nofas8)
saveRDS(all8, "./8_subset_stats_performance.RDS")


