data.reduce <-
  function(data,            # data matrix to reduce based on predictor variable significance, either biological data (type 1, without X and Y columns), environmental or combined (sitepair) dataset
           datatype = "sp", # type of data to reduce: "sp"= site pair; or "pred"= environmental predictor
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

      cdata.sigtest <- rbind(as.matrix(as.logical(c("T","T","T","T","T","T"))),sigtest,sigtest)
      new.data <- data[,cdata.sigtest]
    }

    if(datatype=="pred"){
      cat("Reducing environmental dataset following significance test result\n")
      cat("\n")

      envdata.sigtest <- rbind(as.matrix(as.logical(c("T","T","T"))),sigtest)
      new.data <- data[,envdata.sigtest]
    }

    else{
      stop("Invalid data type!")
    }

    cat("Data reduction finished\n")
    cat("\n")

    return(new.data)
  }
