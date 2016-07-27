#' @title Function to calculate predictor variable drop contributions on GDM model
#' @description ...
#' @param caldata Calibration data
#' @param valdata Validation data
#' @param performance Performance metric to be used ("rmse" or "r2"), set as "r2" per default
#' @param geo Optional use of geographical distance as predictor in GDM model, set as FALSE per default
#' @return Returns model performance value
#' @export

gdm.validate <-
  function(caldata,                # calibration data
           valdata,                # validation data
           performance = "r2",     # performance metric to be used ("rmse" or "r2"), set as "r2" per default
           geo = F)                # optional use of geographical distance as predictor in GDM model set as FALSE per default
  {

    # v.1
    #
    # p. j. leitao - 20th June 2016
    #
    # function to perform independent validation of GDM model
    #
    # delivers model performance value
    #
    # requires gdm
    #

    # dependencies configuration
    #require(gdm)

    # fitting  model

    model <- gdm::gdm(caldata, geo = geo)

    # calculating predicted dissimilarities for partial model

    predicted <- gdm::predict.gdm(model,valdata)
    predicted <- as.data.frame(predicted)

    # calculating model performance
    if (performance == "r2"){
      cat("\n")
      cat("Calculating model R-square (R2)...\n")
      cat("\n")
      performance <- (cor(valdata[,1], predicted)^2)*100
    }

    if (performance == "rmse") {
      cat("\n")
      cat("Calculating model Root Mean Square Error (RMSE)...\n")
      cat("\n")
      performance <- sqrt(mean((valdata[,1]-predicted)^2))
    }

    else{
      stop("Invalid performance type!")
    }

    cat("Model performance calculated\n")
    cat("\n")

    return(performance)
  }

