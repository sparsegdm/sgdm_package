#' @title Function to map NMDS transformed gdm site pair dissimilarity predictions
#' @description Function to predict dissimilarities between sample site pairs using a gdm model. Resulting multidimesional dissimilarities
#'              are transformed to a NMDS (non-metrical multidimensional scaling) feature space with a defined number of components.
#'              If a raster object with the same set of predictor variables as used for gdm model building is provided, the transformed
#'              dissimilarities will be mapped using k-nearest neighbor imputation.
#' @param spData Compiled dataset as output from "format.site.pair".
#' @param predMap Raster object with the same set of predictior variables as used for gdm model building; if no raster object is provided NMDS scores will be returned as data.frame.
#' @param model A gdm.model object for dissimilarity prediction.
#' @param k number of resulting NMDS components (default = 3); if k = 0 number of components will be derived using a mean NMDS stress value threshold (t) after 20 iterations of NMDS using all given predictor variables.
#' @param t NMDS stress value threshold (default = 0.1).
#' @return NMDS object with k number of components based on gdm predicted sample pair dissimilarities if no raster is provided; NMDS raster map with k number of layers.
#' @export



gdm.map<- function(spData,        # site pair table as from formatsitetable function in gdm package
                   predMap,       # Raster object with the same set of predictior variables as used for gdm model building; if no raster object is provided NMDS scores will be returned as data.frame.
                   model,         # gdm.model
                   k = 3,           # number of NMDS components to extract; if not specified number of components will be derived be NMDS stress value
                   t = 0.1)         # NMDS stress value threshold to extract number of components if k is not specified
{
  # v.1
  #
  # m. schwieder; c. senf - 27th July 2016
  #
  # function to map NMDS transformed dissimilarities
  #
  # delivers NMDS outputs as score table or raster map
  #
  # requires gdm, vegan, raster and yaImpute
  #

  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("raster package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("vegan package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("yaImpute", quietly = TRUE)) {
    stop("yaImpute package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # data reading and dependencies configuration
  #require(vegan)
  #require(gdm)
  #require(yaImpute)
  #require(raster)

  pairs <- nrow(spData)
  n2 <- (1+sqrt(1+8*pairs))/2
  nVars <- (ncol(spData)-6)/2
  dummy_ID <- data.frame(ID=1:n2)

  # derive predData from input spData
  first_sample <- cbind(dummy_ID[1,1], spData[1,3:4], spData[1,7:(nVars+6)])
  predData0 <- cbind(dummy_ID[2:n2,1],spData[1:(n2-1),5:6],spData[1:(n2-1),(nVars+7):ncol(spData)])
  names <- c("ID","X","Y", paste0("Pred.", 1:(nVars)))
  colnames(first_sample) <- names
  colnames(predData0) <- names
  predData <- rbind(first_sample, predData0)

  cat("\n")
  cat("Performing gdm mapping\n")
  cat("\n")

  # predict dissimilarities for sample pairs
  sample.pair <- spData
  sample.pair.diss <- gdm::predict.gdm(model, sample.pair)

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
        table_nmds <- vegan::monoMDS(sample.pair.diss.mat, k=i, model ="global", maxit=1000)

        # Stress values of model
        # stressplot(table_nmds)
        stress[j,i]<-table_nmds$stress
      }

      cat("Mean stress value after 20 iterations with k =",i, "is", apply(stress[i],2, mean),"\n")
      cat("\n")
    }

    mean_stress<-as.matrix(apply(stress,2, mean))

    k<- as.numeric(length(which(mean_stress>t))+1)
  }

  # Perform NMDS on sample dissimialirites
  cat("\n")
  cat("Performing NMDS transformation on sample pair sites with",k, "components\n")
  cat("\n")

  sample_nmds <- vegan::monoMDS(sample.pair.diss.mat, k=k, model="global", maxit = 1000)

  cat("\n")
  cat("NMDS stress value:",sample_nmds$stress)
  cat("\n")

  # Extract scores from NMDS model
  nmds_scores <- as.matrix(vegan::scores(sample_nmds))

  # Check if raster stack is provided for map output.
  if (missing(predMap)){
    #return(nmds_scores)
    return(sample_nmds)
  }

  else{
    data.type.check<-class(predMap)

    if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer"){
      cat("\n")
      cat("predMap is data type ", data.type.check, " and will be used to impute the NMDS transformed gdm dissimilarity predictions \n")
      cat("\n")

      image.df <- raster::as.data.frame(predMap, xy=TRUE)
      image.df <- cbind(1:nrow(image.df), image.df)

      # check data consistency between predData and predMap
      if(raster::nlayers(predMap)!= nVars) stop("Raster image must have same numbers of layers as numbers of predictors!")

      names(image.df) <- c("ID", "X", "Y", paste0("Pred.", 1:(nVars)))


      ### Prepare data for imputation

      imputation.df <- cbind(nmds_scores,predData)
      rownames(imputation.df) <- paste0("I", 1:nrow(imputation.df))

      # Create Imputation model
      impute_model <- yaImpute::yai(x = imputation.df[,(k+4):ncol(imputation.df)], y = imputation.df[,1:k], method = "euclidean")

      # apply model to new image pixel
      impute_model_image <- yaImpute::newtargets(impute_model, image.df[, 4:ncol(image.df)])

      # impute image
      impute_image <- yaImpute::impute(impute_model_image)

      # create list of NMDS counts
      r_list <- vector("list", length = k)

      # prepare raster for each NMDS component
      for(i in 1:k){
        r <- raster::subset(predMap, 1)
        raster::values(r) <- impute_image[,i]
        names(r) <- names(impute_image)[i]
        r_list[[i]] <- r
      }

      # stack NMDS raster
      impute_map <- raster::stack(r_list)


    }
  }
  # plot and return NMDS raster
  plot(impute_map, col=gray.colors(100 , start = 0.9, end = 0.1 ))
  return(impute_map)
}

