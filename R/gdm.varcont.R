#' @title Calculates predictor variable (drop) contributions in GDM model
#'
#' @description
#' This function calculates predictor variable (drop) contributions in GDM model.
#'
#' Variable contributions are calculated by the loss in explained deviance in the GDM model, when dropping each individual predictor variable at a time.
#'
#' It requires a combined site pair dataset ("spData" format) and the optional use of geographical distance as predictor variable in the GDM.
#'
#' For more details relating to "spData" data format, check \code{gdm} package.
#'
#' @param spData Combined site pair dataset ("spData" format).
#' @param geo Optional use of geographical distance as predictor in GDM model. Set to \code{FALSE} per default
#' @return Returns variable contribution in percentage
#' @export

gdm.varcont <-
  function(spData,
           geo = FALSE)
  {

    cat("\n")
    cat("Calculating GDM model variable contributions\n")
    cat("\n")

    # data reading

    l1 <- (ncol(spData)-6)/2
    m1a <- l1+6

    varnames <- colnames(spData[7:m1a])

    # creating output object

    contribs <- matrix(0,l1,1)

    # calculating full model

    gdm.mod <- gdm::gdm(spData, geo=geo)

    # iteratively dropping each variable and calculating its drop contribution

    for (h in 1:l1) {
      i1 <- h+6
      i2 <- h+l1+6
      dvdata <- spData[,-i2]
      dvdata <- dvdata[,-i1]
      contribs[h,1]<- (gdm.mod$explained) - ((assign(paste("gdm.drop.",h, sep=""), gdm::gdm(dvdata, geo = geo)))$explained)
    }

    # filling ouput object

    colnames(contribs) <- "contribution"
    rownames(contribs) <- varnames

    cat("GDM model variable contributions calculated\n")
    cat("\n")

    return(contribs)
  }
