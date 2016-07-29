#' @title Reduces predictor data based on significance test
#' @description This function reduces the predictor variables of a dataset based on the significance test, resulting from function gdm.varsig. It can be applied to either predictor or site pair datasets.
#' @param data Data to be reduced based on predictor variable significance. It can be either a predictor dataset ("predData" format) or combined Site pair table ("spData" format).
#' @param datatype Type of data to de reduced: \code{pred} for predictor data ("predData" format) or \code{sp} for site pair data ("spData" format).
#' @param sigtest Predictor variable contribution significance test, as output by gdm.varsig function.
#' @return Returns reduced environmental data in the same format as the input data.
#' @export

data.reduce <-
  function(data,
           datatype = "sp",
           sigtest)
  {

    # Reducing dataset

    if(datatype=="sp"){
      cat("Reducing compiled dataset following significance test result\n")
      cat("\n")

      spData.sigtest0 <- append(as.logical(c("T","T","T","T","T","T")), sigtest)
      spData.sigtest <- append(spData.sigtest0, sigtest)
      new.data <- data[,spData.sigtest]
    }

    if(datatype=="pred"){
      cat("Reducing environmental dataset following significance test result\n")
      cat("\n")

      predData.sigtest <- append(as.logical(c("T","T","T")),sigtest)
      new.data <- data[,predData.sigtest]
    }

    cat("Data reduction finished\n")
    cat("\n")

    return(new.data)
  }
