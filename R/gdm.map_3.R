gdm.map_2 <- function(spData,       # Site pair table as from Formatsitetable gdm function
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
  if (!"yaImpute" %in% installed.packages()){
    stop("Package 'yaImpute' must be installed!")
  }

  # data reading and dependencies configuration
  require(vegan)
  require(gdm)
  require(yaImpute)
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

      image.df <- as.data.frame(predMap, xy=TRUE)
      image.df <- cbind(1:nrow(image.df), image.df) 
      
      if(nlayers(predMap)!= nVars) stop("Raster image must have same numbers of layers as numbers of predictors!")
      
      names(image.df) <- c("ID", "X", "Y", paste0("Pred.", 1:(nVars)))
      
      
      ### Prepare data for imputation
      
      imputation.df <- cbind(nmds_output,predData)
      rownames(imputation.df) <- paste0("I", 1:nrow(imputation.df))
      
      impute_model <- yai(x = imputation.df[,7:ncol(imputation.df)], y = imputation.df[,1:3], method = "euclidean")
      impute_model_image <- newtargets(impute_model, image.df[, 4:ncol(image.df)])
      impute_image <- impute(impute_model_image)
      
      r_list <- vector("list", length = k)
      for(i in 1:k){
        r <- subset(predMap, 1)
        values(r) <- impute_image[,i]
        names(r) <- names(impute_image)[i]
        r_list[[i]] <- r  
      }

      impute_map <- stack(r_list)
    
      return(impute_map)
  }
}
