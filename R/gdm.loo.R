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
