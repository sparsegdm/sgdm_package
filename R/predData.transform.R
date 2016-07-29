#' @title Function to perform SCCA transformation of environmental data
#' @description This function transforms the given environmental input data into sparse canonical components,
#'              based on the sparse canonical vectors as derived from \code {sgdm.best}.Input data can either be provided
#'              as data frame or raster object.Input data must have the same number of predictors as used for deriving the sparse canonical vectors
#' @param predData Environmental data as data frame (Format: Plot_ID, X, Y, env.data) or as a raster object
#' @param v.best sparse canonical vectors as out from \code {sgdm.best}
#' @return Returns CCA transformed environmental data
#' @export

predData.transform <-
  function(predData,
           v.best)
  {

    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("raster package needed for this function to work. Please install it.",
           call. = FALSE)
    }


    data.type.check<-class(predData)

    if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer"){
      cat("\n")
      cat("predMap is data type ", data.type.check, " and will be converted to data frame \n")
      cat("\n")

      # save raster for output
      r <- raster::subset(predData, 1)

      # convert raster
      image.df <- raster::as.data.frame(predData, xy=TRUE)
      image.df <- cbind(1:nrow(image.df), image.df)

      names(image.df) <- c("Plot_ID", "X", "Y", paste0("b", 1:(ncol(image.df)-3)))

      predData<-image.df
    }


    if(nrow(v.best)!= (ncol(predData)-3)) stop("Environmental data must have the same number of predictors as used for deriving sparse canonical vectors!")

    # extracting Plot_ID and coordinates
    gi <- predData[,1:3]

    predData1 <- as.matrix(predData[,4:ncol(predData)])
    predData.new0 <- predData1 %*% v.best
    predData.new <- as.data.frame(cbind(gi,predData.new0))


    if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer"){

      r_list <- vector("list", length = (ncol(predData.new)-3))

      for(i in 1:(ncol(predData.new)-3)){

        raster::values(r) <- predData.new[,i+3]
        names(r) <- names(predData.new)[i+3]
        r_list[[i]] <- r
      }

      # stack NMDS raster
      predData.new<- raster::stack(r_list)
    }

    return(predData.new)

  }


