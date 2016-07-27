#' @title Function to calculate predictor variable drop contributions on GDM model
#' @description ...
#' @param spData Site pair table as from formatsitetable function in gdm package
#' @param geo Optional use of geographical distance as predictor in GDM model, set as FALSE per default
#' @return Returns variable contribution in percentage
#' @export

gdm.varcont <-
  function(spData,              # site pair table as from formatsitetable function in gdm package
           geo = FALSE)         # optional use of geographical distance as predictor in GDM model, set as FALSE per default
  {

    # v.3
    #
    # p. j. leitao - 1st July 2016
    #
    # function to calculate predictor variable drop contributions on GDM model
    #
    # delivers variable contribution in percentage
    #
    # requires gdm
    #

    # data reading and dependencies configuration
    #require(gdm)

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
