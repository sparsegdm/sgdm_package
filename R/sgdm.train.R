#' @title
#' Estimates penalization parameters for SGDM model
#'
#' @description
#' This function estimates the penalization parameters for the SCCA within the SGDM model.
#'
#' The SCCA parameterization is done in a heuristic grid search manner, by testing all possible pairs of parameters (for both the biological and the predictor datasets), according to the resulting GDM performance (RMSE) in 5-fold cross-validation.
#'
#' It requires a predictor dataset ("predData" format), a biological dataset ("bioData" format), the number of components to be extracted in the SCCA, the penalization values to be tested (for both datasets), and the optional use of geographical distance as predictor variable in the GDM.
#'
#' This current implementation only allows biological data in the format 1 using abundance values, as described in the \code{gdm} package.
#'
#' For more details relating to "bioData" and "predData" data formats, check \code{gdm} package.
#'
#' @param predData Predictor dataset ("predData" format).
#' @param bioData Biological dataset ("bioData" format).
#' @param k Number of sparse canonical components to be extracted. Set to 10 per default.
#' @param predPenalization Vector with predictor data penalisation values to be tested in the grid search procedure (between 0 and 1). Set to \code{(0.6, 0.7, 0.8, 0.9, 1)} per default.
#' @param bioPenalization Vector with biological data penalisation values to be tested in the grid search procedure (between 0 and 1). Set to \code{(0.6, 0.7, 0.8, 0.9, 1)} per default.
#' @param geo Optional use of geographical distance as predictor in GDM model. Set to \code{FALSE} per default
#' @return Returns performance matrix with RMSE values for each tested penalization parameter pair.
#' @export

sgdm.train <-
  function(predData,
           bioData,
           k = 10,
           predPenalization = seq(0.6, 1, 0.1),
           bioPenalization = seq(0.6, 1, 0.1),
           geo = F)
    {

    cat("\n")
    cat("Running SGDM model paramerization\n")
    cat("\n")

    # data reading

    j1 <- ncol(predData)

    j2 <- ncol(bioData)
    n2 <- nrow(bioData)
    t2 <- n2-1
    pairc <- n2*t2

    latlong <- as.matrix(predData[,2:3])
    id <- as.matrix(predData[,1])

    r <- as.matrix(bioData[,2:j2])
    p <- as.matrix(predData[,4:j1])

    br <- length(bioPenalization)
    bc <- length(predPenalization)

    # creating output performance matrix
    perf.matrix <- matrix(ncol = bc, nrow = br, data = 0)
    rownames(perf.matrix) = bioPenalization
    colnames(perf.matrix) = predPenalization

    # initialising grid search of SCCA penalization parameters
    cat("Grid search for setting SCCA penalization:\n")

    for (px in bioPenalization) for (pz in predPenalization) {
      cat("\n")
      cat(paste("Penalization on biological data (x) =",px,"; penalization on predcitor data (y) =",pz,"\n"))
      cat("\n")
      cat("SCCA Model:\n")

      # running SCCA
      cca <- PMA::CCA(r, p, typex="standard",typez="standard", penaltyx=px,
                 penaltyz=pz, K=k, niter=50, v=NULL, trace=TRUE, standardize=TRUE,
                 xnames=NULL, znames=NULL)

      # extracting canonical vectors
      v <- cca$v

      # tranforming environmental data into canonical components
      c <- p %*% v
      cgi <- cbind(id,latlong,c)
      cgi <- as.data.frame(cgi)
      colnames(cgi)[1]<- "Plot_ID"

      # compiling dataset
      spData <- gdm::formatsitepair(bioData, 1, dist = "bray", abundance = TRUE,
                              siteColumn = "Plot_ID", XColumn="X",YColumn="Y",
                              predData = cgi)

      # calulating GDM model performance
      result <- gdm.cv(spData,nfolds=5,metric="bray",geo=geo)

      # feeding performance matrix
      perf.matrix[paste(px), paste(pz)] <- result[1]
    }

    cat("\n")
    cat("\n")
    cat("Finished SGDM parameterization: performance raster created\n")
    cat("\n")

    return(perf.matrix)
  }
