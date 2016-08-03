#' @title Performs (independent) validation of GDM model
#'
#' @description
#' This function performs (independent) validation of a GDM model.
#'
#' It requires (independent) calibration and validation datasets in combined site pair ("spData") data format, as well as the definition of the performance measure to be used in the validation and the optional use of geographical distance as predictor variable in the gdm.
#'
#' #' For more details relating to "spData" data format, check \code{gdm} package.
#'
#' @param caldata Calibration dataset, in combined site pair ("spData") data format.
#' @param valdata Validation dataset, in combined site pair ("spData") data format.
#' @param performance Performance metric to be used for validation: \code{"rmse"} for root mean square error (RMSE) or \code{"r2"} for coefficient of determination (r2). Set as \code{"r2"} per default.
#' @param geo Optional use of geographical distance as predictor in GDM model. set to \code{FALSE} per default.
#' @return Returns model performance value.
#' @export

gdm.validate <-
  function(caldata,
           valdata,
           performance = "r2",
           geo = F)
  {

    cat("\n")
    cat("GDM model validation\n")
    cat("\n")

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

