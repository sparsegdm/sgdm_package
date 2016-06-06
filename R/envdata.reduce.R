envdata.reduce <-
  function(envdata,     # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           sigtest)     # predictor variable contribution significance test, as output by gdm.varsig
  {

    #
    # p. j. leit?o - 2nd May 2016
    #
    # function for reducing environmental data based on significance test of GDM predictor variable contributions
    #
    # delivers reduced environmental data matrix
    #

    # Reducing dataset

    cat("Reducing environmental dataset following significance test result\n")
    cat("\n")

    data.sigtest <- rbind(as.matrix(as.logical(c("T","T","T"))),sigtest)
    new.envdata <- envdata[,data.sigtest]

    cat("Data reduction of environmental dataset finished\n")
    cat("\n")

    return(new.envdata)
  }
