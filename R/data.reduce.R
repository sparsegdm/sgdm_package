#' @title Reduces predictor data based on significance test
#'
#' @description
#' This function reduces the predictor variables of a dataset based on the significance test, resulting from function gdm.varsig. It can be applied to either predictor or site pair datasets.
#'
#' For more details relating to "PredData" and "spData" data formats, check \code{gdm} package.
#'
#' @param data Data to be reduced based on predictor variable significance. It can be either a predictor dataset ("predData" format) or combined site pair table ("spData" format).
#' @param datatype Type of data to de reduced: \code{pred} for predictor data ("predData" format) or \code{sp} for site pair data ("spData" format).
#' @param sigtest Predictor variable contribution significance test, as output by \code{gdm.varsig} function.
#'
#' @return Returns reduced environmental or site pair data.
#'
#' @export

data.reduce <-
  function(data,
           datatype = "sp",
           sigtest)
  {

    # Reducing dataset

    if(datatype=="sp"){
      cat("Reducing compiled dataset following significance test result\n")
      cat("\n")

      spData.sigtest0 <- append(as.logical(c("T","T","T","T","T","T")), sigtest)
      spData.sigtest <- append(spData.sigtest0, sigtest)
      new.data <- data[,spData.sigtest]
    }

    if(datatype=="pred"){
      cat("Reducing environmental dataset following significance test result\n")
      cat("\n")


      data.type.check<-class(data)

      if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer"){
        cat("\n")
        cat("Input data is type ", data.type.check,"\n")
        cat("\n")

        # save raster for output
        r <- raster::subset(data, 1)

        # convert raster
        image.df <- raster::as.data.frame(data, xy=TRUE)
        image.df <- cbind(1:nrow(image.df), image.df)

        names(image.df) <- c("Plot_ID", "X", "Y", paste0("b", 1:(ncol(image.df)-3)))

        data<-image.df
      }

      predData.sigtest <- append(as.logical(c("T","T","T")),sigtest)
      new.data <- data[,predData.sigtest]
    }


    if(data.type.check == "RasterStack" || data.type.check == "RasterBrick" || data.type.check == "RasterLayer"){

      r_list <- vector("list", length = (ncol(new.data)-3))

      for(i in 1:(ncol(new.data)-3)){

        raster::values(r) <- new.data[,i+3]
        names(r) <- names(new.data)[i+3]
        r_list[[i]] <- r
      }

      # stack component raster
      new.data<- raster::stack(r_list)
    }

    cat("Data reduction finished\n")
    cat("\n")

    return(new.data)
  }
