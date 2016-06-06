sgdm.best <-
  function(grid.matrix,       # performance matrix as output from sgdm.grid
           envdata,           # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
           biodata,           # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
           output = "g",      # type of output: "g" = gdm model; "c" = sparse canonical components; "v" = sparse canonical vectors; default = gdm model
           comps = 10,        # number of sparce canonical components to be calculated, set as 10 per default
           metric="bray",     # only needed if output = "g"; dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as "bray curtis" per default
           geo = F)           # only needed if output = "g"; optional use of geographical distance as predictor in GDM model, set as FALSE per default

      {

    # v.3
    #
    # p. j. leitao - 6th June 2016
    #
    # function to retrieve the best SGDM model, SCCA canonical components or SCCA canonical vectors, as resulting from the SCCA parameter estimation using sgdm.gridsearch
    #
    # delivers GDM model, sparse canonical components or vectors, extracted with parameter pair with respective best performance
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

    if (output == "v") {
      cat("Sparse canonical vectors created\n")
      cat("\n")

      return (v.best)
    }

    else {
      # transforming environmental data into canonical components

      c.best <- as.matrix(envdata[,4:j1])
      c.best <- c.best %*% v.best
      cgi <- cbind(id,latlong,c.best)
      cgi <- as.data.frame(cgi)
      colnames(cgi)[1] <- "Plot_ID"

      if (output == "c") {
              cat("Sparse canonical components created\n")
      cat("\n")

      return (cgi)
      }

      else {
        # compiling data

        cdata <- data.read(cgi,biodata,metric=metric)

        # running GDM model

        gdm.mod <- gdm(cdata,geo=geo)

        cat("Best SGDM model created\n")
        cat("\n")

        return (gdm.mod)
      }
    }
  }
