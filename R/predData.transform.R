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

    # a3 <- ncol(spData)
    # a0 <- (a3-6)/2
    # a1 <- 6+a0
    # a2 <- a1+1
    # b0 <- nrow(spData)
    j1 <- ncol(predData)
    j2 <- ncol(bioData)

    gi <- predData[,1:3]
    predData0 <- predData[,4:j1]

    # data splitting
    # gi <- spData[,1:6]
    # pair1 <- spData[,7:a1]
    # pair2 <- spData[,a2:a3]

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

    # transforming environmental data into canonical components
    # c.best.pair1 <- as.matrix(pair1)
    # c.best.pair1 <- c.best.pair1 %*% v.best
    # c.best.pair2 <- as.matrix(pair2)
    # c.best.pair2 <- c.best.pair2 %*% v.best

    predData1 <- as.matrix(predData0)
    predData.new0 <- predData1 %*% v.best
    predData.new <- as.data.frame(cbind(gi,predData.new0))

    # spData.new <- cbind(gi,c.best.pair1,c.best.pair2)
    # spData.new <- as.data.frame(spData.new)

    return(predData.new)
  }
