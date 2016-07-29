#' @title
#' Retrieves the best SGDM model, SCCA canonical components or SCCA canonical vectors, as resulting from the SGDM parameter estimation
#'
#' @description
#' This function retrieves the best SGDM model, SCCA canonical components or SCCA canonical vectors, as resulting from the SGDM parameter estimation with the \code{gdm.train} function.
#'
#' The parameter pair with the lowest RMSE value is selected to run the SCCA on the biological and predictor datasets. If \code{output} = "m" delivers the GDM model built on the extracted SCCA components; if \code{output} = "c" delivers the SCCA components that result in the best GDM model; and if If \code{output} = "v" delivers the SCCA canonical vectors used to tranform the predictor data into the canonical components.
#'
#' It requires a performance matrix as resulting from the \code{gdm.train} function, a predictor dataset ("predData" format), a biological dataset ("bioData" format), the type of output, the number of components to be extracted in the SCCA and the optional use of geographical distance as predictor variable in the GDM.
#'
#' This current implementation only allows biological data in the format 1 using abundance values, as described in the \code{gdm} package.
#'
#' For more details relating to "bioData" and "predData" data formats, check \code{gdm} package.
#'
#' @param perf.matrix Performance matrix as output from \code{sgdm.train} function.
#' @param predData Predictor dataset ("predData" format).
#' @param bioData Biological dataset ("bioData" format).
#' @param output Type of output: "m" = gdm model; "c" = sparse canonical components; "v" = sparse canonical vectors; Set as "m" per default.
#' @param k Number of sparce canonical components to be calculated, set as 10 per default
#' @param geo only needed if output = "m"; optional use of geographical distance as predictor in GDM model, set as FALSE per default
#' @return Returns a GDM model, the sparse canonical components, or the sparse canonical vectors, depending on the output defined. The default is \code{output} = "m", which returns a GDM model object.
#' @export

sgdm.best <-
  function(perf.matrix,
           predData,
           bioData,
           output = "m",
           k = 10,
           geo = F)
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

    cca.best <- PMA::CCA(bioData[,2:j2], predData[,4:j1], typex="standard",typez="standard", penaltyx=cname,
                    penaltyz=rname, K=k, niter=1000, v=NULL, trace=TRUE, standardize=TRUE,
                    xnames=NULL, znames=NULL)

    # extracting canonical vectors

    v.best <- cca.best$v

    if (output == "v") {

      cat("\n")
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

        cat("\n")
        cat("Sparse canonical components created\n")
        cat("\n")

      return (cgi)
      }

      if (output == "m") {
        # compiling data

        spData <- gdm::formatsitepair(bioData, 1, dist = "bray", abundance = TRUE,
                                 siteColumn = "Plot_ID", XColumn = "X",YColumn = "Y",
                                 predData = cgi)

        # running GDM model

        gdm.mod <- gdm::gdm(spData,geo=geo)

        cat("\n")
        cat("Best SGDM model created\n")
        cat("\n")

        return (gdm.mod)
      }
    }
  }
