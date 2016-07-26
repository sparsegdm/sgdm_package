gdm.map_2 <- function(spData,        # Site pair table as from Formatsitetable gdm function
                    predMap,       # Data frame that has the same predictior variables as used for gdm model building
                    model,         # gdm.model
                    k=3,           # number of NMDS components to extract; if not specified number of components will be derived be NMDS stress value
                    t=0.1)         # NMDS stress value threshold to extract number of components if k is not specified
{

  # checking dependencies
  if (!"gdm" %in% installed.packages()){
    stop("Package 'gdm' must be installed!")
  }
  if (!"vegan" %in% installed.packages()){
    stop("Package 'vegan' must be installed!")
  }
  if (!"raster" %in% installed.packages()){
    stop("Package 'raster' must be installed!")
  }

  # data reading and dependencies configuration
  require(vegan)
  require(gdm)
  require(raster)

  pairs <- nrow(spData)
  n2 <- (1+sqrt(1+8*pairs))/2
  nVars<-(ncol(spData)-6)/2
  dummy_ID<-data.frame(ID=1:n2)

  # derive predData from input spData
  first_sample<-cbind(dummy_ID[1,1], spData[1,3:4], spData[1,7:(nVars+6)])
  predData0 <- cbind(dummy_ID[2:n2,1],spData[1:(n2-1),5:6],spData[1:(n2-1),(nVars+7):ncol(spData)])
  names<-c("ID","X","Y", paste0("Pred.", 1:(nVars)))
  colnames(first_sample)<-names
  colnames(predData0)<-names
  predData<-rbind(first_sample, predData0)

  cat("\n")
  cat("Performing gdm mapping\n")
  cat("\n")

  # predict dissimilarities for sample pairs
  sample.pair <- spData
  sample.pair.diss<-predict.gdm(model, sample.pair)

  #sample.pair.diss.mat<-as.data.frame(sample.pair.diss)

  X <- diag(x=0, nrow(predData),nrow(predData))
  X[upper.tri(X, diag=FALSE)] <- sample.pair.diss
  X <- X + t(X) - diag(diag(X))

  sample.pair.diss.mat<-as.matrix(X)

  if (k == 0){
    cat("\n")
    cat("Derive number of NMDS components based on stress values below" ,t, "\n")
    cat("\n")

    # derive k for NMDS based on 20 iterations and all possible components
    # stress: > 0.05 excellent
    #         > 0.1  great
    #         > 0.2  good/ok
    #         > 0.3  poor

    stress<-matrix(nrow=20,ncol=nrow(sample.pair.diss.mat))
    stress<-as.data.frame(stress)

    for (i in 1:nrow(sample.pair.diss.mat)){
      for (j in 1:20){
        table_nmds <- monoMDS(sample.pair.diss.mat, k=i, model ="global", maxit=1000)

        # Stress values of model
        # stressplot(table_nmds)
        stress[j,i]<-table_nmds$stress
      }

      cat("Mean stress value after 20 iterations with k =",i, "is", apply(stress[i],2, mean),"\n")
      cat("\n")
    }

    mean_stress<-as.matrix(apply(stress,2, mean))
    # plot(mean_stress, main="Mean NMDS stress values out of 20 iterations")
    k<- as.numeric(length(which(mean_stress>t))+1)
  }

  # Perform NMDS on sample dissimialirites
  cat("\n")
  cat("Performing NMDS transformation on sample pair sites with",k, "components\n")
  cat("\n")

  sample_nmds <- monoMDS(sample.pair.diss.mat, k=k, model="global", maxit = 1000)

  cat("\n")
  cat("NMDS stress value:",sample_nmds$stress)
  cat("\n")

  # stressplot(sample_nmds)


  # Extract scores from NMDS model
  nmds_scores <- as.matrix(scores(sample_nmds))

  # calculate nmds components for predicted sample dissimilarities
        # sample.sample<-matrix(sample.pair.diss, ncol=nrow(predData), byrow=FALSE)
        # sample_nmds_comp<-sample.pair.diss.mat %*% nmds_scores

  # Check if raster stack is provided for map output.
  if (missing(predMap)){
    return(nmds_scores)
  }

  else{
    data.type.check<-class(predMap)

    if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer" ){
      cat("\n")
      cat("predMap is a raster object and will be transformed to a matrix\n")
      cat("\n")

      rst<-predMap

      image.df <- as.data.frame(predMap)
      rst.ID <- matrix(nrow = nrow(image.df), ncol = 1, data = c(1:nrow(image.df)))

      # create image subset and transform to spatialPointsDF to derive XY coordinates
      # img.sub<-subset(rst,1)
      # coords<-rasterToPoints(img.sub)
      # names(img.sub.df)<-c("x","y","1")
      # x<-coords[,1]
      # y<-coords[,2]

      # using dummy xy coordinates
      x <- matrix(nrow = nrow(image.df), ncol = 1, data = 0)
      y <- matrix(nrow = nrow(image.df), ncol = 1, data = 0)

      image<-as.data.frame(cbind(rst.ID,x,y,image.df))
      names.mt <- as.list(c("Plot_ID","X","Y", paste0("MapPred.", 1:(ncol(image.df)))))
      colnames(image) <- names.mt
      predMap<-as.data.frame(image)
    }

            # prepare sample/map site pair table
            # check if both input data fit each other

            # cat("\n")
            # cat("Checking consistency between map and sample data\n")
            # cat("\n")
            #
            # j1 <- ncol(predData)
            # j3 <- ncol(predMap)
            # n.sites <- nrow(predData)
            # weight<-1
            #
            # if (j1 != j3){
            #   stop("Map and samples not consistent, check input data!\n")
            # }
            #
            # cat("\n")
            # cat("Data check OK! Compiling dataset...\n")
            # cat("\n")
            #
            # # creat site pair table
            # map <- vector(mode = "list", length = n.sites)
            # map <- lapply(1:length(map), function(x)map[[x]] <- predMap)
            # map <- do.call("rbind", map)
            #
            # sites <- vector(mode = "list", length = n.sites)
            # sites <- lapply(1:n.sites, function(x)do.call("rbind", rep(list(predData[x,]), nrow(predMap))))
            # sites <- do.call("rbind", sites)
            #
            # map.sample.pair <- cbind(rep(1, nrow(map)), rep(weight, nrow(map)), map[,2:3], sites[,2:3], map[,4:ncol(map)], sites[,4:ncol(sites)])
            # names(map.sample.pair) <- c("distance", "weight", "p.x", "p.y", "s.x", "s.y", paste0("p.env.", 1:(ncol(map)-3)), paste0("s.env.", 1:(ncol(sites)-3)))

    # create biodata dummy with as many rows as pixel in the prediction map
    dummy<-as.data.frame(matrix(data=c(1:10), ncol=3, nrow=nrow(predMap)))
    colnames(dummy)<-c("Plot_ID","sp1","sp2")
    dummy[,1] <- predMap[,1]
    dummy[,2] <- round(runif(nrow(predMap), 1, 10))
    dummy[,3] <- round(runif(nrow(predMap), 1, 10))

    # prepare map/map site pair table
    map.pair<-formatsitepair(dummy,1,siteColumn = "Plot_ID", XColumn="X",YColumn="Y",predData = predMap)

    # predict dissimilarities for map sample site table based on gdm model
    map.predict <- predict(gdm.model,map.pair)

    # Prepare dissimilarity vector as dist object for NMDS transformation
    X <- diag(x=0, nrow(image.df),nrow(image.df))
    X[upper.tri(X, diag=FALSE)] <- map.predict
    X <- X + t(X) - diag(diag(X))
    map.predict.mt<-as.matrix(X)


    # Perform NMDS on sample dissimialirites
    cat("\n")
    cat("Performing NMDS transformation on map pair sites with",k, "components\n")
    cat("\n")

    # Transform predictions using NMDS
    map.predict.nmds <- monoMDS(map.predict.mt, k=k, model="global", maxit = 1000)


    cat("\n")
    cat("NMDS stress value:",map.predict.nmds$stress)
    cat("\n")


    # Extract NMDS scores
    map_trans<-scores(map.predict.nmds)

   # names(map_trans) <- c(paste0("NMDS_",1:(ncol(map_trans))))

    if(data.type.check == "matrix" || data.type.check == "data.frame") {
      cat("\n")
      cat("Done with dissimilarity prediction and NMDS transformation. Output type R dataFrame.")
      cat("\n")

      return(map_trans)
    }

    # derive raster extent from input and assign values to single raster layers

    for (i in 1:ncol(map_trans)){
      # create base raster and assign the first NMDS component values
      if(i == 1){
        nmds_map<-raster(rst, layer=1)
        values(nmds_map)<-map_trans[,1]
        name_1<-colnames(map_trans)
        names(nmds_map)<-name_1[1]
      }

      else{
        # add raster layers with following NMDS values assigned
        nmds_rst<- paste("nmds_rst", i, sep = "")
        assign(nmds_rst, raster(rst, layer=1))

        nmds<- paste("nmds", i, sep = "")
        assign(nmds, as.matrix(map_trans[,i]))
        nmds_values<-get(nmds)
        rst_temp<-get(nmds_rst)
        values(rst_temp)<-nmds_values
        names(rst_temp)<-name_1[i]

        # stack NMDS component raster layers
        nmds_map<-stack(nmds_map,rst_temp)

        plot(nmds_map)
      }
    }

    if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer" ){
      cat("\n")
      cat("Finished dissimilarity prediction and NMDS transformation. Output type R RasterStack\n")
      cat("\n")

      return(nmds_map)
    }
  }
}
