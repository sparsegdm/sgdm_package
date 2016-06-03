# data.read  <-
#   function(envdata,            # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
#            biodata,            # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
#            metric="bray",      # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence)
#            abundance = TRUE)    # weighting of sample plots, set as 1 per default. when given must be a vector as long as the sample size
#   {
#
#     # v.2
#     #
#     # p. j. leit?o - 2nd May 2016
#     #
#     # function to compile species and environmental data matrices into single dataset, suitable for input in GDM4tables
#     #
#     # delivers compiled dataset
#     #
#     # requires gdm
#     #
#
#     # data reading and consistency check
#
#     require(gdm)
#
#     cat("\n")
#     cat("Checking data consistency: ")
#
#     j1 <- ncol(envdata)
#     k1 <- j1-1
#     l1 <- k1-2
#     m1a <- l1+6
#     m1b <- m1a+1
#     m1c <- m1a+l1
#     n1 <- nrow(envdata)
#     ndcols <- 6+(l1*2)
#     varnames <- colnames(envdata[4:j1])
#
#     j2 <- ncol(biodata)
#     k2 <- j2-1
#     n2 <- nrow(biodata)
#     pairs <- ((n2^2)-n2)/2
#
#     if (n1 == 1){
#       cat("ERROR! Data has one sample only, no pairwise combinations possible! ")
#       cat("Please note that an empty object may have been created!\n")
#       cat("\n")
#     }
#
#     else{
#       if (n1 != n2){
#         cat("ERROR! Biological and environmental data not consistent, check input data!\n")
#         cat("\n")}
#
#       else{
#         cat("Data check OK!\n")
#         cat("\n")
#
#         cat(paste("Dataset with",k2,"taxa,",l1,"environmental variables and",pairs,"pairwise combinations"))
#         cat("\n")
#         cat("\n")
#
#         cat(paste("Distance metric used:",metric))
#         cat("\n")
#         cat("\n")
#
#         cdata <- formatsitepair(biodata, 1, dist = metric, abundance = abundance,
#                                 siteColumn = "Plot_ID", XColumn="X",YColumn="Y",
#                                 predData = envdata)
#
#         cat("Dataset loaded")
#         cat("\n")
#         cat("\n")
#
#         return(cdata)
#       }
#     }
#   }


# gdm.loo <-
#   function(cdata,                  # compiled dataset as output from "data.read" function
#            biodata,                # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
#            metric="bray",          # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence)
#            performance="rmse",     # performance metric to be used ("rmse" or "r2"), set as "rmse" per default
#            geo = F)                # optional use of geographical distance as predictor in GDM model set as FALSE per default
#   {
#
#     # v.2
#     #
#     # p. j. leit?o - 2nd May 2014
#     #
#     # function to perform leave-one-out cross-validation of GDM model
#     #
#     # delivers model performance value
#     #
#     # requires gdm, vegan
#     #
#
#     # data reading and dependencies configuration
#
#     require(gdm)
#     require(vegan)
#
#     j2 <- ncol(biodata)
#     k2 <- j2-1
#     n2 <- nrow(biodata)
#     t2 <- n2-1
#     u2 <- (2*t2)-1
#     pairc <- n2*t2
#     pairs <- (n2*t2)/2
#
#     #   consistency check
#
#     cat("\n")
#     cat("GDM Leave-one-out cross-validation\n")
#     cat("\n")
#
#     cat("Checking data consistency: ")
#
#     if (pairs != nrow(cdata)){
#       cat("ERROR! Inconsistent number of samples! ")
#       cat("\n")
#     }
#
#     else{
#       cat("Data check OK!\n")
#       cat("\n")
#
#       # creating performance test data
#
#       r <- biodata[,2:j2]
#
#       perf.test <- matrix (ncol = 2, nrow = pairc, data = 0)
#       colnames(perf.test) = c("observed", "predicted")
#
#       # selecting sample to hold out
#
#       a <- 1                #1 at start
#       s2 <- t2              #n-1 at start
#
#       for (h in 1:n2) {
#         cat(paste("Held-out sample",h,"of",n2,"\n"))
#
#         index <- as.data.frame(matrix(ncol = n2, nrow = n2, data = 1))
#         index[h,] <- 0
#         index[,h] <- 0
#         indexd <- as.dist(index)
#
#         index.sel <- as.data.frame(matrix(ncol = 1, nrow = pairs, data = 0))
#
#         o <- 1
#         q <- 1
#         for (p in 2:n2){
#           for (i in p:n2){
#             index.sel[q,1] <- indexd[q]
#             q <- q+1
#           }
#           o <- o+1
#         }
#
#         # defining calibration and validation datasets
#
#         caldata <- cdata[which(index.sel==1),]
#         valdata <- cdata[which(index.sel==0),]
#
#         r.sample1 <- as.data.frame(matrix(ncol = k2, nrow = t2, data = 0))
#
#         for (i in 1:t2) {
#           r.sample1[i,1:k2] <- as.data.frame(r[h,])
#         }
#
#         r.sample2 <- r[-h,]
#
#         colnames(r.sample1) <- colnames(r.sample2)
#         r.samples <- rbind(r.sample1,r.sample2)
#
#         # calculating observed dissimilarities
#
#         observed <- vegdist(r.samples, method=metric)
#         observed <- observed[t2:u2]
#         observed <- as.data.frame(observed)
#
#         # running partial model
#
#         partial.gdm <- gdm(caldata, geo = geo)
#
#         # calculating predicted dissimilarities for partial model
#
#         predicted <- predict(partial.gdm,valdata)
#         predicted <- as.data.frame(predicted)
#
#         perf.test[a:s2,1] <- observed[1:t2,1]
#         perf.test[a:s2,2] <- predicted[1:t2,1]
#
#         a <- a+t2
#         s2 <- s2+t2
#       }
#
#       # calculating model performance
#
#       if (performance == "r2"){
#         cat("\n")
#         cat("Calculating model R-square (R2)...\n")
#         cat("\n")
#         performance <- (cor(perf.test[,1], perf.test[,2])^2)*100
#       }
#
#       if (performance == "rmse") {
#         cat("\n")
#         cat("Calculating model Root Mean Square Error (RMSE)...\n")
#         cat("\n")
#         performance <- sqrt(mean((perf.test[,1]-perf.test[,2])^2))
#       }
#
#       cat("Model performance calculated\n")
#       cat("\n")
#
#       return (performance)
#     }
#   }


# gdm.varcont <-
#   function(cdata,               # compiled dataset as output from "data.read" function
#            geo = FALSE)         # optional use of geographical distance as predictor in GDM model, set as FALSE per default
#   {
#
#     # v.2
#     #
#     # p. j. leit?o - 2nd May 2016
#     #
#     # function to calculate predictor variable drop contributions on GDM model
#     #
#     # delivers variable contribution in percentage
#     #
#     # requires gdm
#     #
#
#     # data reading
#
#     cat("Calculating GDM model variable contributions\n")
#     cat("\n")
#
#     l1 <- (ncol(cdata)-6)/2
#     m1a <- l1+6
#
#     varnames <- colnames(cdata[7:m1a])
#
#     # creating output object
#
#     contribs <- matrix(0,l1,1)
#
#     # calculating full model
#
#     gdm.mod <- gdm(cdata, geo=geo)
#
#     # iteratively dropping each variable and calculating its drop contribution
#
#     for (h in 1:l1) {
#       i1 <- h+6
#       i2 <- h+l1+6
#       dvdata <- cdata[,-i2]
#       dvdata <- dvdata[,-i1]
#       contribs[h,1]<- (gdm.mod$explained) - ((assign(paste("gdm.drop.",h, sep=""), gdm(dvdata, geo = geo)))$explained)
#     }
#
#     # filling ouput object
#
#     colnames(contribs) <- "contribution"
#     rownames(contribs) <- varnames
#
#     cat("GDM model variable contributions calculated\n")
#     cat("\n")
#
#     return(contribs)
#   }


gdm.varsig <-
  function(envdata,         # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           biodata,         # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
           metric = "bray", # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as "bray curtis" per default
           geo = F,         # optional use of geographical distance as predictor in GDM model, set as FALSE per default
           perm = 100,      # number of matrix permutations, set as 100 per default
           sig = 0.05)      # significance level, set as 0.05 per default
  {

    # v.2
    #
    # p. j. leit?o - 2nd October 2016
    #
    # function to test statistical significance predictor variable drop contributions on GDM model
    #
    # requires gdm, vegan, ecodist
    #
    # delivers results of significance test for each predictor variable
    #

    # dependencies

    require(vegan)

    # data reading

    j2 <- ncol(biodata)
    biodata10000 <- cbind(biodata[,1],(trunc(biodata[,2:j2]*10000)))
    colnames(biodata10000)[1]<-colnames(biodata)[1]

    cdata <- data.read(envdata,biodata10000,metric=metric)

    cat("Testing for significance of GDM model variable contributions\n")
    cat("\n")

    l1 <- (ncol(cdata)-6)/2

    j2 <- ncol(biodata)
    n2 <- nrow(biodata)

    w <- perm
    w1 <- w+1

    # creating variable contrveibution matrix

    contribs <- matrix(0,l1,w1)

    # initiating significance test: species matrix permutations

    cat("Checking statistical significance of GDM variable contributions:\n")
    cat("Performing species matrix permutations\n")
    cat("\n")

    bio.m <- as.matrix(biodata10000[,2:j2])
    bio.perm <- permatfull(bio.m, fixedmar = "row", times = w)

    perm.m <- matrix(0,w,1)

    # calculating variable contributions of original model

    cat("Running original GDM model and calculating variable contributions\n")
    cat("\n")

    contribs[,1] <- gdm.varcont(cdata, geo = geo)

    # calculating variable contributions of permuted models

    cat("Calculating variable contributions of permuted models\n")
    cat("\n")

    for (t in 1:w) {
      cat(paste("Permutation",t,"of",w,"\n", sep = " "))
      cat("\n")

      q <- t+1
      dist <- vegdist(as.data.frame(bio.perm$perm[t]), method=metric)
      pdata <- cdata
      pdata[,1] <- dist

      contribs[,q] <- gdm.varcont(pdata, geo = geo)
    }

    # Feeding variable contribution matrix

    cont.perm <- matrix(0,l1,w)

    for (g in 1:l1){
      for (f0 in 1:w){
        f <- f0+1
        cont.perm[g,f0] <- contribs[g,1] < contribs[g,f]
      }
    }

    cont.sig <- matrix(0,l1,1)

    # significance test

    for (e in 1:l1){
      cont.sig[e,] <- sum(cont.perm[e,])/w
    }

    cat("Significance test of GDM model variable contributions completed\n")
    cat("\n")

    return(cont.sig<=sig)
  }



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


cdata.reduce <-
  function(cdata,               # compiled dataset as output from "data.read" function
           sigtest)             # predictor variable contribution significance test, as output by gdm.varsig
  {

    #
    # p. j. leit?o - 12th May 2014
    #
    # function for reducing compiled dataset based on significance test of GDM predictor variable contributions
    #
    # delivers reduced compiled dataset
    #

    # Reducing dataset

    cat("Reducing compiled dataset following significance test result\n")
    cat("\n")

    cdata.sigtest <- rbind(as.matrix(as.logical(c("T","T","T","T","T","T"))),sigtest,sigtest)
    new.cdata <- cdata[,cdata.sigtest]

    cat("Data reduction of compiled dataset finished\n")
    cat("\n")

    return(new.cdata)
  }


sgdm.gridsearch <-
  function(envdata,                   # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           biodata,                   # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
           comps = 10,                # number of sparce canonical components to be calculated, set as 10 per default
           metric="bray",             # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as bray curtis" per default
           env.penalization=penalize, # vector with possible penalisation values to be applied on the environmental data matrix (between 0 and 1)
           bio.penalization=penalize, # vector with possible penalisation values to be applied on the biological data matrix (between 0 and 1)
           geo = F)                   # optional use of geographical distance as predictor in GDM model, set as FALSE per default
  {

    # v.2
    #
    # p. j. leit?o - 2nd May 2016
    #
    # function to perform parameter estimation of SCCA via grid search, based on GDM leave one out cross-validated performances (RMSE)
    #
    # delivers performance matrix with RMSE values for each SCCA penalization parameter pair
    #

    # dependencies

    require(PMA)

    # data reading

    cat("\n")
    cat("Running SGDM model paramerization\n")
    cat("\n")

    j1 <- ncol(envdata)

    j2 <- ncol(biodata)
    n2 <- nrow(biodata)
    t2 <- n2-1
    pairc <- n2*t2

    latlong <- as.matrix(envdata[,2:3])
    id <- as.matrix(envdata[,1])

    r <- as.matrix(biodata[,2:j2])
    p <- as.matrix(envdata[,4:j1])

    penalize <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

    br <- length(bio.penalization)
    bc <- length(env.penalization)

    # creating output performance matrix

    perfmatrix <- matrix(ncol = bc, nrow = br, data = 0)
    rownames(perfmatrix) = bio.penalization
    colnames(perfmatrix) = env.penalization

    # creating performance test

    perf.test <- matrix (ncol = 2, nrow = pairc, data = 0)
    colnames(perf.test) = c("observed", "predicted")

    # initialising grid search of SCCA penalization parameters

    cat("Grid search for setting SCCA penalization:\n")

    for (px in bio.penalization) for (pz in env.penalization) {
      cat("\n")
      cat(paste("Penalization on species data (x) =",px,"; penalization on environmental data (y) =",pz,"\n"))
      cat("\n")
      cat("SCCA Model:\n")

      # running SCCA

      cca <- CCA(r, p, typex="standard",typez="standard", penaltyx=px,
                 penaltyz=pz, K=comps, niter=50, v=NULL, trace=TRUE, standardize=TRUE,
                 xnames=NULL, znames=NULL)

      # extracting canonical vectors

      v <- cca$v

      # tranforming environmental data into canonical components

      c <- p %*% v
      cgi <- cbind(id,latlong,c)
      cgi <- as.data.frame(cgi)
      colnames(cgi)[1]<- "Plot_ID"

      # compiling dataset

      cdata <- data.read(cgi,biodata,metric=metric)

      # calulating GDM model performance

      performance <- gdm.loo(cdata,biodata,metric=metric,geo=geo)

      # feeding performance matrix

      perfmatrix[paste(px), paste(pz)] <- performance
    }

    cat("\n")
    cat("\n")
    cat("Finished SGDM parameterization: performance raster created\n")
    cat("\n")

    return(perfmatrix)
  }



sgdm.best <-
  function(grid.matrix,           # performance matrix as output from sgdm.grid
           envdata,               # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           biodata,               # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
           comps = 10,            # number of sparce canonical components to be calculated, set as 10 per default
           metric="bray",         # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as "bray curtis" per default
           geo = F)               # optional use of geographical distance as predictor in GDM model, set as FALSE per default
  {

    # v.2
    #
    # p. j. leit?o - 2nd May 2016
    #
    # function to retrieve the best SGDM model, as resulting from the SCCA parameter estimation using sgdm.grid
    #
    # delivers GDM model using sparce canonical components extracted with parameter pair with respective best performance
    #

    # dependencies

    require(PMA)

    # data reading

    cat("Retrieving best SGDM model after parameterization\n")
    cat("\n")

    j1 <- ncol(envdata)
    j2 <- ncol(biodata)

    latlong <- as.matrix(envdata[,2:3])
    id <- as.matrix(envdata[,1])

    # reading SCCA parameterization from performance matrix

    min.index<-arrayInd(which.min(grid.matrix),dim(grid.matrix))
    rname <- as.numeric(rownames(grid.matrix)[min.index[1]])
    cname <- as.numeric(colnames(grid.matrix)[min.index[2]])

    # running SCCA

    cca.best <- CCA(biodata[,2:j2], envdata[,4:j1], typex="standard",typez="standard", penaltyx=cname,
                    penaltyz=rname, K=comps, niter=1000, v=NULL, trace=TRUE, standardize=TRUE,
                    xnames=NULL, znames=NULL)

    # extracting canonical vectors

    v.best <- cca.best$v

    # transforming environmental data into canonical components

    c.best <- as.matrix(envdata[,4:j1])
    c.best <- c.best %*% v.best
    cgi <- cbind(id,latlong,c.best)
    cgi <- as.data.frame(cgi)
    colnames(cgi)[1]<- "Plot_ID"

    # compiling data

    cdata <- data.read(cgi,biodata,metric=metric)

    # running GDM model

    gdm.mod <- gdm(cdata,geo=geo)

    cat("Best SGDM model created\n")
    cat("\n")

    return (gdm.mod)
  }


scc.best <-
  function(grid.matrix,       # performance matrix as output from sgdm.grid
           envdata,           # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           biodata,           # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
           comps = 10)        # number of sparce canonical components to be calculated, set as 10 per default
  {

    #
    # p. j. leit?o - 9th May 2014
    #
    # function to retrieve the SCCA canonical components corresponding to the best SGDM model, as resulting from the SCCA parameter estimation using sgdm.grid
    #
    # delivers sparse canonical components extracted with parameter pair with respective best performance
    #

    # dependencies

    require(PMA)

    # data reading

    cat("Retrieving sparse canonical components corresponding to the best SGDM model after parameterization\n")
    cat("\n")

    j1 <- ncol(envdata)
    j2 <- ncol(biodata)

    latlong <- as.matrix(envdata[,2:3])
    id <- as.matrix(envdata[,1])

    # reading SCCA parameterization from performance matrix

    min.index<-arrayInd(which.min(grid.matrix),dim(grid.matrix))
    rname <- as.numeric(rownames(grid.matrix)[min.index[1]])
    cname <- as.numeric(colnames(grid.matrix)[min.index[2]])

    # running SCCA

    cca.best <- CCA(biodata[,2:j2], envdata[,4:j1], typex="standard",typez="standard", penaltyx=cname,
                    penaltyz=rname, K=comps, niter=1000, v=NULL, trace=TRUE, standardize=TRUE,
                    xnames=NULL, znames=NULL)

    # extracting canonical vectors

    v.best <- cca.best$v

    # transforming environmental data into canonical components

    c.best <- as.matrix(envdata[,4:j1])
    c.best <- c.best %*% v.best
    cgi <- cbind(id,latlong,c.best)
    cgi <- as.data.frame(cgi)
    names(cgi)[1] <- "Plot_ID"

    cat("Sparse canonical components created\n")
    cat("\n")

    return (cgi)
  }


scv.best <-
  function(grid.matrix,       # performance matrix as output from sgdm.grid
           envdata,           # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           biodata,           # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
           comps = 10)        # number of sparce canonical components to be calculated, set as 10 per default
  {

    #
    # p. j. leit?o - 9th May 2014
    #
    # function to retrieve the SCCA canonical vectors corresponding to the best SGDM model, as resulting from the SCCA parameter estimation using sgdm.grid
    #
    # delivers sparse canonical vectors extracted with parameter pair with respective best performance
    #

    # dependencies

    require(PMA)

    # data reading

    cat("Retrieving sparse canonical vectors corresponding to the best SGDM model after parameterization\n")
    cat("\n")

    j1 <- ncol(envdata)
    j2 <- ncol(biodata)

    latlong <- as.matrix(envdata[,2:3])
    id <- as.matrix(envdata[,1])

    # reading SCCA parameterization from performance matrix

    min.index<-arrayInd(which.min(grid.matrix),dim(grid.matrix))
    rname <- as.numeric(rownames(grid.matrix)[min.index[1]])
    cname <- as.numeric(colnames(grid.matrix)[min.index[2]])

    # running SCCA

    cca.best <- CCA(biodata[,2:j2], envdata[,4:j1], typex="standard",typez="standard", penaltyx=cname,
                    penaltyz=rname, K=comps, niter=1000, v=NULL, trace=TRUE, standardize=TRUE,
                    xnames=NULL, znames=NULL)

    # extracting canonical vectors

    v.best <- cca.best$v

    cat("Sparse canonical vectors created\n")
    cat("\n")

    return (v.best)
  }
