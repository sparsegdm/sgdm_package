gdm.varcont <-
  function(cdata,               # compiled dataset as output from "data.read" function
           geo = FALSE)         # optional use of geographical distance as predictor in GDM model, set as FALSE per default
  {

    # v.2
    #
    # p. j. leit?o - 2nd May 2016
    #
    # function to calculate predictor variable drop contributions on GDM model
    #
    # delivers variable contribution in percentage
    #
    # requires gdm
    #

    # data reading

    cat("Calculating GDM model variable contributions\n")
    cat("\n")

    l1 <- (ncol(cdata)-6)/2
    m1a <- l1+6

    varnames <- colnames(cdata[7:m1a])

    # creating output object

    contribs <- matrix(0,l1,1)

    # calculating full model

    gdm.mod <- gdm(cdata, geo=geo)

    # iteratively dropping each variable and calculating its drop contribution

    for (h in 1:l1) {
      i1 <- h+6
      i2 <- h+l1+6
      dvdata <- cdata[,-i2]
      dvdata <- dvdata[,-i1]
      contribs[h,1]<- (gdm.mod$explained) - ((assign(paste("gdm.drop.",h, sep=""), gdm(dvdata, geo = geo)))$explained)
    }

    # filling ouput object

    colnames(contribs) <- "contribution"
    rownames(contribs) <- varnames

    cat("GDM model variable contributions calculated\n")
    cat("\n")

    return(contribs)
  }
