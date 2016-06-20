# data.read  <-
#   function(envdata,            # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
#            biodata,            # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
#            metric="bray",      # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence)
#            abundance = TRUE)    # weighting of sample plots, set as 1 per default. when given must be a vector as long as the sample size
#   {
#
#     # v.2
#     #
#     # p. j. leit?o - 2nd May 2016
#     #
#     # function to compile species and environmental data matrices into single dataset, suitable for input in GDM4tables
#     #
#     # delivers compiled dataset
#     #
#     # requires gdm
#     #
#
#     # data reading and consistency check
#
#     require(gdm)
#
#     cat("\n")
#     cat("Checking data consistency: ")
#
#     j1 <- ncol(envdata)
#     k1 <- j1-1
#     l1 <- k1-2
#     m1a <- l1+6
#     m1b <- m1a+1
#     m1c <- m1a+l1
#     n1 <- nrow(envdata)
#     ndcols <- 6+(l1*2)
#     varnames <- colnames(envdata[4:j1])
#
#     j2 <- ncol(biodata)
#     k2 <- j2-1
#     n2 <- nrow(biodata)
#     pairs <- ((n2^2)-n2)/2
#
#     if (n1 == 1){
#       stop("Data has one sample only, no pairwise combinations possible!")
#     }
#
#     else{
#       if (n1 != n2){
#         stop("Biological and environmental data not consistent, check input data!")
#       }
#
#       else{
#         cat("Data check OK!\n")
#         cat("\n")
#
#         cat(paste("Dataset with",k2,"taxa,",l1,"environmental variables and",pairs,"pairwise combinations"))
#         cat("\n")
#         cat("\n")
#
#         cat(paste("Distance metric used:",metric))
#         cat("\n")
#         cat("\n")
#
#         cdata <- formatsitepair(biodata, 1, dist = metric, abundance = abundance,
#                                 siteColumn = "Plot_ID", XColumn="X",YColumn="Y",
#                                 predData = envdata)
#
#         cat("Dataset loaded")
#         cat("\n")
#         cat("\n")
#
#         return(cdata)
#       }
#     }
#   }
