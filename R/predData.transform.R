predData.transform <-
  function(predData,
           bioData,
           perf.matrix,
           comps = 10)
  {

    # checking dependencies
    if (!"PMA" %in% installed.packages()){
      stop("Package 'PMA' must be installed!")
    }

    # data reading and dependencies configuration
    require(PMA)

    j1 <- ncol(predData)
    j2 <- ncol(bioData)

    gi <- predData[,1:3]
    predData0 <- predData[,4:j1]

    # reading SCCA parameterization from performance matrix
    cat("Retrieving best SGDM model after parameterization\n")
    cat("\n")

    min.index<-arrayInd(which.min(perf.matrix),dim(perf.matrix))
    rname <- as.numeric(rownames(perf.matrix)[min.index[1]])
    cname <- as.numeric(colnames(perf.matrix)[min.index[2]])

    # running SCCA
    cca.best <- CCA(bioData[,2:j2], predData0, typex="standard",typez="standard", penaltyx=cname,
                    penaltyz=rname, K=comps, niter=1000, v=NULL, trace=TRUE, standardize=TRUE,
                    xnames=NULL, znames=NULL)

    # extracting canonical vectors
    v.best <- cca.best$v

    predData1 <- as.matrix(predData0)
    predData.new0 <- predData1 %*% v.best
    predData.new <- as.data.frame(cbind(gi,predData.new0))

    return(predData.new)
  }
