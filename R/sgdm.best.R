#' Function to retrieve the best SGDM model, SCCA canonical components or SCCA canonical vectors, as resulting from the SCCA parameter estimation using sgdm.gridsearch
#'
#' @param perf.matrix Performance matrix as output from sgdm.grid
#' @param predData Environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
#' @param bioData Biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
#' @param output Type of output: "m" = gdm model; "c" = sparse canonical components; "v" = sparse canonical vectors; default = gdm model
#' @param comps Number of sparce canonical components to be calculated, set as 10 per default
#' @param metric Only needed if output = "m"; dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as "bray curtis" per default
#' @param geo only needed if output = "m"; optional use of geographical distance as predictor in GDM model, set as FALSE per default
#' @return Returns GDM model, sparse canonical components, or sparse canonical vectors (default is \code{output = "m"}, which returns a GDM model object)
#' @export

sgdm.best <-
  function(perf.matrix,       # performance matrix as output from sgdm.grid
           predData,          # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           bioData,           # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
           output = "m",      # type of output: "m" = gdm model; "c" = sparse canonical components; "v" = sparse canonical vectors; default = gdm model
           comps = 10,        # number of sparce canonical components to be calculated, set as 10 per default
           metric="bray",     # only needed if output = "m"; dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as "bray curtis" per default
           geo = F)           # only needed if output = "m"; optional use of geographical distance as predictor in GDM model, set as FALSE per default

      {

    # v.3
    #
    # p. j. leitao - 6th June 2016
    #
    # function to retrieve the best SGDM model, SCCA canonical components or SCCA canonical vectors, as resulting from the SCCA parameter estimation using sgdm.gridsearch
    #
    # delivers GDM model, sparse canonical components or vectors, extracted with parameter pair with respective best performance
    #

    # data reading and dependencies configuration
    #require(gdm)
    #require(PMA)

    cat("Retrieving sparse canonical components corresponding to the best SGDM model after parameterization\n")
    cat("\n")

    j1 <- ncol(predData)
    j2 <- ncol(bioData)

    latlong <- as.matrix(predData[,2:3])
    id <- as.matrix(predData[,1])

    # reading SCCA parameterization from performance matrix

    min.index<-arrayInd(which.min(perf.matrix),dim(perf.matrix))
    rname <- as.numeric(rownames(perf.matrix)[min.index[1]])
    cname <- as.numeric(colnames(perf.matrix)[min.index[2]])

    # running SCCA

    cca.best <- CCA(bioData[,2:j2], predData[,4:j1], typex="standard",typez="standard", penaltyx=cname,
                    penaltyz=rname, K=comps, niter=1000, v=NULL, trace=TRUE, standardize=TRUE,
                    xnames=NULL, znames=NULL)

    # extracting canonical vectors

    v.best <- cca.best$v

    if (output == "v") {
      cat("Sparse canonical vectors created\n")
      cat("\n")

      return (v.best)
    }

    else {
      # transforming environmental data into canonical components

      c.best <- as.matrix(predData[,4:j1])
      c.best <- c.best %*% v.best
      cgi <- cbind(id,latlong,c.best)
      cgi <- as.data.frame(cgi)
      colnames(cgi)[1] <- "Plot_ID"

      if (output == "c") {
              cat("Sparse canonical components created\n")
      cat("\n")

      return (cgi)
      }

      if (output == "m") {
        # compiling data

        spData <- formatsitepair(bioData, 1, dist = metric, abundance = TRUE,
                                 siteColumn = "Plot_ID", XColumn = "X",YColumn = "Y",
                                 predData = cgi)

        # running GDM model

        gdm.mod <- gdm(spData,geo=geo)

        cat("Best SGDM model created\n")
        cat("\n")

        return (gdm.mod)
      }
    }
  }
