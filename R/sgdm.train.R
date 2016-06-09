sgdm.train <-
  function(envdata,                   # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           biodata,                   # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
           comps = 10,                # number of sparce canonical components to be calculated, set as 10 per default
           metric="bray",             # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as bray curtis" per default
           env.penalization=penalize, # vector with possible penalisation values to be applied on the environmental data matrix (between 0 and 1)
           bio.penalization=penalize, # vector with possible penalisation values to be applied on the biological data matrix (between 0 and 1)
           geo = F)                   # optional use of geographical distance as predictor in GDM model, set as FALSE per default
  {

    # v.2
    #
    # p. j. leit?o - 2nd May 2016
    #
    # function to perform parameter estimation of SCCA via grid search, based on GDM leave one out cross-validated performances (RMSE)
    #
    # delivers performance matrix with RMSE values for each SCCA penalization parameter pair
    #

    # dependencies

    require(PMA)

    # data reading

    cat("\n")
    cat("Running SGDM model paramerization\n")
    cat("\n")

    j1 <- ncol(envdata)

    j2 <- ncol(biodata)
    n2 <- nrow(biodata)
    t2 <- n2-1
    pairc <- n2*t2

    latlong <- as.matrix(envdata[,2:3])
    id <- as.matrix(envdata[,1])

    r <- as.matrix(biodata[,2:j2])
    p <- as.matrix(envdata[,4:j1])

    penalize <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

    br <- length(bio.penalization)
    bc <- length(env.penalization)

    # creating output performance matrix

    perfmatrix <- matrix(ncol = bc, nrow = br, data = 0)
    rownames(perfmatrix) = bio.penalization
    colnames(perfmatrix) = env.penalization

    # creating performance test

    perf.test <- matrix (ncol = 2, nrow = pairc, data = 0)
    colnames(perf.test) = c("observed", "predicted")

    # initialising grid search of SCCA penalization parameters

    cat("Grid search for setting SCCA penalization:\n")

    for (px in bio.penalization) for (pz in env.penalization) {
      cat("\n")
      cat(paste("Penalization on species data (x) =",px,"; penalization on environmental data (y) =",pz,"\n"))
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

      cdata <- data.read(cgi,biodata,metric=metric)

      # calulating GDM model performance

      performance <- gdm.cv(cdata,nfolds=5,metric=metric,geo=geo)

      # feeding performance matrix

      perfmatrix[paste(px), paste(pz)] <- performance
    }

    cat("\n")
    cat("\n")
    cat("Finished SGDM parameterization: performance raster created\n")
    cat("\n")

    return(perfmatrix)
  }
