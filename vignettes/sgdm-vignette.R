## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  devtools::install_git(...)
#  library(sgdm)

## ---- echo=FALSE---------------------------------------------------------
library(sgdm)

## ---- hide=TRUE----------------------------------------------------------
sgdm.gs <- sgdm.train(predData = spectra, bioData = trees, k = 30, predPenalization = seq(0.6, 1, 0.1), bioPenalization = seq(0.6, 1, 0.1), geo = F)

## ------------------------------------------------------------------------
print(sgdm.gs, digits = 4)

## ---- hide=TRUE----------------------------------------------------------
sgdm.sccbest <- sgdm.best(perf.matrix = sgdm.gs, predData = spectra, bioData = trees, output = "c", k = 30)

## ------------------------------------------------------------------------
spData.best <- gdm::formatsitepair(trees, 1, dist = "bray", abundance = TRUE,
               siteColumn = "Plot_ID", XColumn = "X",YColumn = "Y",
               predData = sgdm.sccbest)

## ------------------------------------------------------------------------
sgdm.model <- gdm::gdm(spData.best)

## ------------------------------------------------------------------------
sgdm.model$explained

## ------------------------------------------------------------------------
plot(sgdm.model$predicted, sgdm.model$observed, xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1)

## ---- hide=TRUE----------------------------------------------------------
sigtest.sgdm <- gdm.varsig(predData = sgdm.sccbest, bioData = trees)

## ------------------------------------------------------------------------
spData.red.best <- data.reduce(spData.best, datatype = "sp", sigtest.sgdm)
sgdm.red.model <- gdm::gdm(spData.red.best)

sgdm.red.model$explained

plot(sgdm.red.model$predicted,sgdm.red.model$observed,xlim=c(0,1),ylim=c(0,1))
abline(0,1)

## ------------------------------------------------------------------------
gdm.cv(spData.best, nfolds = 10, performance = "rmse")

## ------------------------------------------------------------------------
gdm.cv(spData.best, nfolds = 10, performance = "r2")

## ------------------------------------------------------------------------
raster::plotRGB(spectral.image, r = 20, g = 10, b = 2, stretch = "hist")

## ------------------------------------------------------------------------
v <- sgdm.best(perf.matrix = sgdm.gs, predData = spectra, bioData = trees, output = "v", k = 30)
component.image <- predData.transform(predData = spectral.image, v = v)
component.image.red <- data.reduce(component.image, datatype = "pred", sigtest = sigtest.sgdm)

## ------------------------------------------------------------------------
mapsgdm.red <- gdm.map(spData.red.best, component.image.red, sgdm.red.model, k = 3)
raster::plotRGB(mapsgdm.red, r = 3, g = 2, b = 1, stretch = "hist")

