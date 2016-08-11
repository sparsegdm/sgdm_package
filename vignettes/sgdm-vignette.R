## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  devtools::install_git("sparsegdm/sgdm_package")
#  library(sgdm)

## ---- echo=FALSE---------------------------------------------------------
library(sgdm)

## ---- results='hide'-----------------------------------------------------
sgdm.gs <- sgdm.param(predData = spectra, bioData = trees, k = 30, predPenalization = seq(0.6, 1, 0.1), bioPenalization = seq(0.6, 1, 0.1), geo = F)

## ------------------------------------------------------------------------
print(sgdm.gs, digits = 4)

## ---- results='hide'-----------------------------------------------------
sgdm.model <- sgdm.best(perf.matrix = sgdm.gs, predData = spectra, bioData = trees, output = "m", k = 30)

## ------------------------------------------------------------------------
summary(sgdm.model)

sgdm.model$explained

## ---- results='hide'-----------------------------------------------------
sgdm.sccbest <- sgdm.best(perf.matrix = sgdm.gs, predData = spectra, bioData = trees, output = "c", k = 30)
sgdm.vbest <- sgdm.best(perf.matrix = sgdm.gs, predData = spectra, bioData = trees, output = "v", k = 30)

## ---- results='hide'-----------------------------------------------------
sigtest.sgdm <- gdm.varsig(predData = sgdm.sccbest, bioData = trees)

## ------------------------------------------------------------------------
sgdm.sccbest.red <- data.reduce(data = sgdm.sccbest, datatype = "pred", sigtest = sigtest.sgdm)

## ------------------------------------------------------------------------
spData.sccabest.red <- gdm::formatsitepair(bioData = trees, bioFormat = 1, dist = "bray",
                                           abundance = TRUE, siteColumn = "Plot_ID", 
                                           XColumn = "X",YColumn = "Y", predData = sgdm.sccbest.red)

sgdm.model.red <- gdm::gdm(data = spData.sccabest.red)

## ------------------------------------------------------------------------
summary(sgdm.model.red)

sgdm.model.red$explained

plot(sgdm.model.red$predicted, sgdm.model.red$observed, xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1)

## ------------------------------------------------------------------------
gdm.cv(spData = spData.sccabest.red, nfolds = 10)

## ------------------------------------------------------------------------
gdm.cv(spData = spData.sccabest.red, nfolds = 10, performance = "r2")

## ------------------------------------------------------------------------
community.samples <- gdm.map(spData = spData.sccabest.red, model = sgdm.model.red, k = 0, t = 0.1)

## ------------------------------------------------------------------------
raster::plotRGB(spectral.image, r = 43, g = 22, b = 12, stretch = "hist")

## ------------------------------------------------------------------------
component.image <- predData.transform(predData = spectral.image, v = sgdm.vbest)
component.image.red <- data.reduce(component.image, datatype = "pred", sigtest = sigtest.sgdm)

## ------------------------------------------------------------------------
map.sgdm.red <- gdm.map(spData = spData.sccabest.red, predMap = component.image.red, model = sgdm.model.red, k = 3)
raster::plotRGB(map.sgdm.red, r = 3, g = 2, b = 1, stretch = "hist")

