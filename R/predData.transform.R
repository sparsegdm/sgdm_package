#' @title Performs SCCA transformation of environmental data
#'
#' @description
#' This function transforms a predictor dataset into sparse canonical components, based on the sparse canonical vectors extracted by function \code{sgdm.best} using output = \code{"v"}.
#'
#' The input data can either be provided as dataframe or raster object, which must have the same number of predictors as used for deriving the sparse canonical vectors. The output data will be delivered in the same format (dataframe or raster object) as the input.
#'
#' If predictor dataset is a raster object, this function requires the previous instalation of the \code{raster} package.
#'
#' For more details relating to "predData" and "spData" data formats, check \code{gdm} package.
#'
#' @param predData Predictor dataset (predData format) or as a raster object.
#' @param v Sparse canonical vectors as extracted by \code{sgdm.best} using output = \code{"v"}.
#' @return Returns environmental data transformed into sparse canonical components for further use with GDM. This dataset is delivered in the same format (data frame or raster object) as the input data.
#' @export

predData.transform <-
  function(predData,
           v)
  {

    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("raster package needed for this function to work. Please install it.",
           call. = FALSE)
    }

    data.type.check<-class(predData)

    if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer"){
      cat("\n")
      cat("predData is data type ", data.type.check,"\n")
      cat("\n")

      # save raster for output
      r <- raster::subset(predData, 1)

      # convert raster
      image.df <- raster::as.data.frame(predData, xy=TRUE)
      image.df <- cbind(1:nrow(image.df), image.df)

      names(image.df) <- c("Plot_ID", "X", "Y", paste0("b", 1:(ncol(image.df)-3)))

      predData<-image.df
    }


    if(nrow(v)!= (ncol(predData)-3)) stop("Environmental data must have the same number of predictors as used for deriving sparse canonical vectors!")

    # extracting Plot_ID and coordinates
    gi <- predData[,1:3]

    predData1 <- as.matrix(predData[,4:ncol(predData)])
    predData.new0 <- predData1 %*% v
    predData.new <- as.data.frame(cbind(gi,predData.new0))


    if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer"){

      r_list <- vector("list", length = (ncol(predData.new)-3))

      for(i in 1:(ncol(predData.new)-3)){

        raster::values(r) <- predData.new[,i+3]
        names(r) <- names(predData.new)[i+3]
        r_list[[i]] <- r
      }

      # stack component raster
      predData.new<- raster::stack(r_list)
    }

    return(predData.new)

  }


