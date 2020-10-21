library(tidyr)
library(dplyr)
library(ggplot2)

plottingDir <- "../../data/texture_results/plots/"

model <- c("Pix", "Pix_FAS", "Pix_FAS_HOS", "FAS", "HOS")

######## BSD data ##########
bsdFile <- "../../data/BSD_results/BSD_Coefs.Rds"
bsdData <- readRDS(bsdFile) %>%
  dplyr::mutate(., Pix = as.integer(Pix), FAS = as.integer(FAS),
                HOS = as.integer(HOS)) %>%
  dplyr::select(., model, everything()) %>%
  as_tibble(.)

# extract parameter names
parNames <- names(bsdData)
indPix <- c(grep("pix", parNames), grep("^LP", parNames))
namesPix <- parNames[indPix]
indFAS <- c(grep("mm", parNames), grep("acr", parNames))
namesFAS <- parNames[indFAS]
indsHOS <- which(!(parNames %in% c(namesPix, namesFAS)))
indsHOS <- indsHOS[-c(1:2)]
namesHOS <- parNames[indsHOS]

###############

### Pixel params ###
pixData <- dplyr::filter(bsdData, Pix == 1) %>%
  dplyr::select(., c(1:5), all_of(namesPix)) %>%
  tidyr::pivot_longer(., namesPix, names_to = "parameter", values_to = "coef")

#### plot results when using pixel data ####
pixelCoefPlot <- pixData %>%
  ggplot(aes(x = model, y = abs(coef), color = model)) +
  geom_point() +
  facet_wrap(~ parameter, ncol = round(sqrt(length(namesPix))), scales = "free") +
  theme_bw()


#### FAS params ######
fasData <- dplyr::filter(bsdData, FAS == 1) %>%
  dplyr::select(., c(1:5), all_of(namesFAS)) %>%
  tidyr::pivot_longer(., namesFAS, names_to = "parameter", values_to = "coef")

cutPoint <- round(length(namesFAS)/3)
namesFAS_split[[1]] <- namesFAS[1:cutPoint]
namesFAS_split[[2]] <- namesFAS[(cutPoint+1):(cutPoint*2)]
namesFAS_split[[3]] <- namesFAS[cutPoint*2+1:length(namesFAS)]

#### plot results when using pixel data ####
fasPlot <- fasData %>%
  ggplot(aes(x = model, y = coef, color = model)) +
  geom_point() +
  facet_wrap(~ parameter, ncol = round(sqrt(length(namesFAS))), scales = "free") +
  theme_bw()






