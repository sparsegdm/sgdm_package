#' Function for reducing environmental data based on significance test of GDM predictor variable contributions
#'
#' @param data Data matrix to reduce based on predictor variable significance, either biological data (type 1, without X and Y columns), environmental or combined (sitepair) dataset
#' @param datatype Type of data to reduce: \code{sp} for site pair or \code{pred} for environmental predictors
#' @param sigtest Predictor variable contribution significance test, as output by \code{\link{gdm.varsig}}
#' @return Returns reduced environmental data matrix
#' @examples
#' data.reduce(...)
#' @export

data.reduce <-
  function(data,            # data matrix to reduce based on predictor variable significance, either biological data (type 1, without X and Y columns), environmental or combined (sitepair) dataset
           datatype = "sp", # type of data to reduce: "sp"= site pair; or "pred"= environmental predictors
           sigtest)         # predictor variable contribution significance test, as output by gdm.varsig
  {

    #
    # p. j. leitao - 2nd May 2016
    #
    # function for reducing environmental data based on significance test of GDM predictor variable contributions
    #
    # delivers reduced environmental data matrix
    #

    # Reducing dataset

    if(datatype=="sp"){
      cat("Reducing compiled dataset following significance test result\n")
      cat("\n")

      # spData.sigtest <- rbind(as.matrix(as.logical(c("T","T","T","T","T","T"))),sigtest,sigtest)
      spData.sigtest0 <- append(as.logical(c("T","T","T","T","T","T")),sigtest.sgdm)
      spData.sigtest <- append(spData.sigtest0, sigtest.sgdm)
      new.data <- data[,spData.sigtest]
    }

    if(datatype=="pred"){
      cat("Reducing environmental dataset following significance test result\n")
      cat("\n")

      # predData.sigtest <- rbind(as.matrix(as.logical(c("T","T","T"))),sigtest)
      predData.sigtest <- append(as.logical(c("T","T","T")),sigtest.sgdm)
      new.data <- data[,predData.sigtest]
    }

    # else{
    #   stop("Invalid data type!")
    # }

    cat("Data reduction finished\n")
    cat("\n")

    return(new.data)
  }
