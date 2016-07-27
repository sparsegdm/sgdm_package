#' Function to perform n-fold cross-validation of gdm model
#'
#' @param spData Compiled dataset as output from "format.site.pair" function
#' @param nfolds Number of folds for cross-validation; default = 10; if nfolds = samples then performs leave-one-out cross-validation
#' @return ABC
#' @export

gdm.cv <-
  function(spData,                 # compiled dataset as output from "data.read" function
           nfolds=10,              # number of folds for cross-validation; default = 10; if nfolds = samples then performs leave-one-out cross-validation
           metric="bray",          # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence)
           performance="rmse",     # performance metric to be used ("rmse" or "r2"), set as "rmse" per default
           geo = F)                # optional use of geographical distance as predictor in GDM model set as FALSE per default
  {
    # v.1
    #
    # p. j. leitao - 6th May 2016
    #
    # function to perform n-fold cross-validation of GDM model
    #
    # delivers model performance value
    #
    # requires gdm
    #

    # data reading and dependencies configuration
    #require(gdm)

    pairs <- nrow(spData)
    n2 <- (1+sqrt(1+8*pairs))/2
    t2 <- n2-1
    pairc <- n2*t2

    cat("\n")
    cat("GDM cross-validation\n")
    cat("\n")

    if (nfolds > n2){
      cat("\n")
      cat("ERROR: Incorrect number of folds for GDM cross-validation\n")
      cat("\n")
    }

    else{

      if (nfolds == n2){
        cat("Performing leave-one-out cross-validation")
        cat("\n")
        cat("\n")

        # creating performance test data
        perf.test <- matrix (ncol = 2, nrow = pairc, data = 0)
        colnames(perf.test) = c("observed", "predicted")

        # selecting sample to hold out
        a <- 1                #1 at start
        s2 <- t2              #n-1 at start

        for (h in 1:n2) {
          cat(paste("Held-out sample",h,"of",n2,"\n"))

          index <- as.data.frame(matrix(ncol = n2, nrow = n2, data = 1))
          index[h,] <- 0
          index[,h] <- 0
          indexd <- as.dist(index)
          index.sel <- as.data.frame(matrix(ncol = 1, nrow = pairs, data = 0))

          o <- 1
          q <- 1
          for (p in 2:n2){
            for (i in p:n2){
              index.sel[q,1] <- indexd[q]
              q <- q+1
            }
            o <- o+1
          }

          # defining calibration and validation datasets
          caldata <- spData[which(index.sel==1),]
          valdata <- spData[which(index.sel==0),]

          observed <- as.data.frame(valdata[,1])

          # running partial model
          partial.gdm <- gdm(caldata, geo = geo)

          # calculating predicted dissimilarities for partial model
          predicted <- predict(partial.gdm,valdata)
          predicted <- as.data.frame(predicted)

          perf.test[a:s2,1] <- observed[1:t2,1]
          perf.test[a:s2,2] <- predicted[1:t2,1]

          a <- a+t2
          s2 <- s2+t2
        }
      }

      if (nfolds < n2){
        cat(paste("Performing ",nfolds,"-fold cross-validation",sep=""))
        cat("\n")
        cat("\n")

        # creating performance test data
        perf.test <- matrix (ncol = 2, nrow = nfolds, data = 0)
        colnames(perf.test) = c("observed", "predicted")

        # selecting fold to hold out
        selector <- rep(seq(1, nfolds, by = 1), length = n2)
        selector <- selector[order(runif(n2, 1, 100))]

        a <- 1                #1 at start

        for (h in 1:nfolds) {
          cat(paste("Fold",h,"of",nfolds,"\n"))

          index <- as.data.frame(matrix(ncol = n2, nrow = n2, data = 1))
          index[selector == h,] <- 0
          index[,selector == h] <- 0
          indexd <- as.dist(index)

          index.sel <- as.data.frame(matrix(ncol = 1, nrow = pairs, data = 0))

          o <- 1
          q <- 1
          for (p in 2:n2){
            for (i in p:n2){
              index.sel[q,1] <- indexd[q]
              q <- q+1
            }
            o <- o+1
          }

          # defining calibration and validation datasets
          caldata <- spData[which(index.sel==1),]
          valdata <- spData[which(index.sel==0),]

          observed <- as.data.frame(valdata[,1])

          # running partial model
          partial.gdm <- gdm(caldata, geo = geo)

          # calculating predicted dissimilarities for partial model
          predicted <- predict(partial.gdm,valdata)
          predicted <- as.data.frame(predicted)

          t3 <- nrow(valdata)
          s2 <- t3

          if (h == 1){
            perf.test <- cbind(observed,predicted)
            colnames(perf.test) = c("observed", "predicted")
          }

          else {
            temp <- cbind(observed,predicted)
            colnames(temp)<-colnames(perf.test)
            perf.test <- rbind(perf.test,temp)
          }

          a <- a+t3
          s2 <- s2+t3
        }
      }

      # calculating model performance
      if (performance == "r2"){
        cat("\n")
        cat("Calculating model R-square (R2)...\n")
        cat("\n")
        performance <- (cor(perf.test[,1], perf.test[,2])^2)*100
      }

      if (performance == "rmse") {
        cat("\n")
        cat("Calculating model Root Mean Square Error (RMSE)...\n")
        cat("\n")
        performance <- sqrt(mean((perf.test[,1]-perf.test[,2])^2))
      }

      cat("Model performance calculated\n")
      cat("\n")

      return (performance)
    }
  }
