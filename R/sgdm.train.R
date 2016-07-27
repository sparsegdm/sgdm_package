#' Function to perform parameter estimation of SCCA via grid search, based on GDM leave one out cross-validated performances (RMSE)
#'
#' @param predData Environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
#' @param bioData Biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
#' @param comps Number of sparce canonical components to be calculated, set as 10 per default
#' @param metric Dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as bray curtis" per default
#' @param predPenalization Vector with possible penalisation values to be applied on the environmental data matrix (between 0 and 1)
#' @param bioPenalization Vector with possible penalisation values to be applied on the biological data matrix (between 0 and 1)
#' @param geo Optional use of geographical distance as predictor in GDM model, set as FALSE per default
#' @return Returns performance matrix with RMSE values for each SCCA penalization parameter pair
#' @examples
#' sgdm.train(...)
#' @export

sgdm.train <-
  function(predData,                          # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           bioData,                           # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
           comps = 10,                        # number of sparce canonical components to be calculated, set as 10 per default
           metric = "bray",                   # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as bray curtis" per default
           predPenalization = seq(0.6, 1, 0.1), # vector with possible penalisation values to be applied on the environmental data matrix (between 0 and 1)
           bioPenalization = seq(0.6, 1, 0.1), # vector with possible penalisation values to be applied on the biological data matrix (between 0 and 1)
           geo = F)                           # optional use of geographical distance as predictor in GDM model, set as FALSE per default
  {

    # v.2
    #
    # p. j. leitao - 2nd May 2016
    #
    # function to perform parameter estimation of SCCA via grid search, based on GDM leave one out cross-validated performances (RMSE)
    #
    # delivers performance matrix with RMSE values for each SCCA penalization parameter pair
    #

    # checking dependencies
    if (!"gdm" %in% installed.packages()){
      stop("Package 'gdm' must be installed!")
    }
    if (!"PMA" %in% installed.packages()){
      stop("Package 'PMA' must be installed!")
    }

    # data reading and dependencies configuration
    require(gdm)
    require(PMA)

    cat("\n")
    cat("Running SGDM model paramerization\n")
    cat("\n")

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
      cca <- CCA(r, p, typex="standard",typez="standard", penaltyx=px,
                 penaltyz=pz, K=comps, niter=50, v=NULL, trace=TRUE, standardize=TRUE,
                 xnames=NULL, znames=NULL)

      # extracting canonical vectors
      v <- cca$v

      # tranforming environmental data into canonical components
      c <- p %*% v
      cgi <- cbind(id,latlong,c)
      cgi <- as.data.frame(cgi)
      colnames(cgi)[1]<- "Plot_ID"

      # compiling dataset
      spData <- formatsitepair(bioData, 1, dist = metric, abundance = TRUE,
                              siteColumn = "Plot_ID", XColumn="X",YColumn="Y",
                              predData = cgi)

      # calulating GDM model performance
      result <- gdm.cv(spData,nfolds=5,metric=metric,geo=geo)

      # feeding performance matrix
      perf.matrix[paste(px), paste(pz)] <- result[1]
    }

    cat("\n")
    cat("\n")
    cat("Finished SGDM parameterization: performance raster created\n")
    cat("\n")

    return(perf.matrix)
  }
