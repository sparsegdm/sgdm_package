cdata.transform <-
  function(cdata,
           endata,
           spdata,
           grid.matrix,
           comps = 10)
  {


    # data reading

    a3 <- ncol(cdata)
    a0 <- (a3-6)/2
    a1 <- 6+a0
    a2 <- a1+1
    b0 <- nrow(cdata)
    j1 <- ncol(endata)
    j2 <- ncol(spdata)

    # data splitting

    gi <- cdata[,1:6]
    pair1 <- cdata[,7:a1]
    pair2 <- cdata[,a2:a3]

    # reading SCCA parameterization from performance matrix

    cat("Retrieving best SGDM model after parameterization\n")
    cat("\n")

    min.index<-arrayInd(which.min(grid.matrix),dim(grid.matrix))
    rname <- as.numeric(rownames(grid.matrix)[min.index[1]])
    cname <- as.numeric(colnames(grid.matrix)[min.index[2]])

    # running SCCA

    cca.best <- CCA(spdata[,2:j2], endata[,4:j1], typex="standard",typez="standard", penaltyx=cname,
                    penaltyz=rname, K=comps, niter=1000, v=NULL, trace=TRUE, standardize=TRUE,
                    xnames=NULL, znames=NULL)

    # extracting canonical vectors

    v.best <- cca.best$v

    # transforming environmental data into canonical components

    c.best.pair1 <- as.matrix(pair1)
    c.best.pair1 <- c.best.pair1 %*% v.best
    c.best.pair2 <- as.matrix(pair2)
    c.best.pair2 <- c.best.pair2 %*% v.best

    cdata.new <- cbind(gi,c.best.pair1,c.best.pair2)
    cdata.new <- as.data.frame(cdata.new)

    return(cdata.new)
  }
