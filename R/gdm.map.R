#' @title
#' Maps community composition as a result of non-linear transformation of the predicted GDM dissimilarities
#'
#' @description
#' This function maps community composition as a result of the non-linear transformation of the predicted GDM dissimilarities.
#'
#' A GDM model is run on the input site pair dataset ("spData" format), in order to derive the predicted dissimilarities between all possible pairs of sites.
#'
#' A non-metrical multidimensional scaling (NMDS) is then applied to transform the predicted dissimilarities to extract the main patterns of community composition within the sampled sites. The number of NMDS components to be extracted can be mannually defined by the user (set to 3 per default). Alternatively, if the number of desired components given is set to 0 (zero), the resulting number of components will be automatically defined based on the resulting stress values (of the NMDS) according to a set threshold \code{t} (set to 0.1 per default).
#'
#' It requires a It requires a combined site pair dataset ("spData" format), a GDM model object and the number of components (\code{k}) to be extracted in the NMDS. If \code{k} is set to 0 then requires the stress value threshold \code{t} for automatically defining the number of NMDS components to be extracted.
#'
#' If a map with predictor variables is given (either as raster object or dataframe), these will be assigned (mapped) along the NMDS axes, by comparing the map with the sampled sites (site pair data), following to a k-nearest neighbor imputation approach. The output map will be delivered in the same format as the input map. For this method to work, the input map must have exactly the same predictor variables as those used for the GDM model, and those in the site pair dataset. This option requires the previous instalation of the \code{yalmpute} package, and in case a raster object is given as input map it also requires the \code{raster} package.
#'
#' For more details relating to "spData" data format, check \code{gdm} package.
#'
#' @param spData Combined site pair dataset ("spData" format). This dataset must include the same predictor variables as those used to build the GDM model below.
#' @param predMap Map with the same prediction variables as used to compile the site pair dataset. This map may be a raster object or a dataframe. The output map will be delivered in the same format as the input predMap. If no variable is given here, the function will only map the community composition on the sampled sites. i.e. the NMDS scores for the site pair data.
#' @param model A GDM model object for dissimilarity prediction.
#' @param k Number of NMDS components to be extracted. Set to 3 per default. If \code{k} = 0 then the number of components to be extracted will be automatically defined following the mean NMDS stress value threshold \code{t} after 20 iterations.
#' @param t NMDS stress value threshold. Set to 0.1 per default. This variable is only valid if \code{k} (above) is set to 0 (zero).
#' @return If a predMap is given returns a NMDS map with k number of layers. If no predMap is given returns NMDS scores for k number of components based on the GDM predicted dissimilarities for the site pair data.
#' @export



gdm.map<- function(spData,
                   predMap,
                   model,
                   k = 3,
                   t = 0.1)
  {

  if (!requireNamespace("raster", quietly = TRUE)) {
    stop("raster package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("yaImpute", quietly = TRUE)) {
    stop("yaImpute package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # data reading

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
    return(nmds_scores)
    # return(sample_nmds)
  }

  else{
    data.type.check<-class(predMap)

    if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer"){
      cat("\n")
      cat("predMap, is data type", data.type.check, "and will be used to impute the NMDS-transformed GDM predicted dissimilarities\n")
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

