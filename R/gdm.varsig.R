# gdm.varsig <-
#   function(envdata,         # environmental data matrix, with Plot_ID, X, Y as three first columns followed by predictor values per plot
#            biodata,         # biological data matrix, with Plot_ID as first column, followed by species occurrence / abundance per plot
#            metric = "bray", # dissimilarity metric to be used ("bray curtis" for abundance or "Jaccard" for presence-absence), set as "bray curtis" per default
#            geo = F,         # optional use of geographical distance as predictor in GDM model, set as FALSE per default
#            perm = 100,      # number of matrix permutations, set as 100 per default
#            sig = 0.05)      # significance level, set as 0.05 per default
#   {
#
#     # v.3
#     #
#     # p. j. leitao - 5th june 2016
#     #
#     # function to test statistical significance predictor variable drop contributions on GDM model
#     #
#     # requires gdm, vegan, ecodist
#     #
#     # delivers results of significance test for each predictor variable
#     #
#
#     # dependencies
#
#     require(vegan)
#
#     # data reading
#
#     j2 <- ncol(biodata)
#
#     cdata <- data.read(envdata,biodata,metric=metric)
#
#     cat("Testing for significance of GDM model variable contributions\n")
#     cat("\n")
#
#     l1 <- (ncol(cdata)-6)/2
#
#     j2 <- ncol(biodata)
#     n2 <- nrow(biodata)
#
#     w <- perm
#     w1 <- w+1
#
#     # creating variable contrveibution matrix
#
#     contribs <- matrix(0,l1,w1)
#
#     # initiating significance test: species matrix permutations
#
#     cat("Checking statistical significance of GDM variable contributions:\n")
#     cat("Performing species matrix permutations\n")
#     cat("\n")
#
#     bio.m <- as.matrix(biodata[,2:j2])
#     bio.perm <- permatfull(bio.m, fixedmar = "row", times = w)
#
#     perm.m <- matrix(0,w,1)
#
#     # calculating variable contributions of original model
#
#     cat("Running original GDM model and calculating variable contributions\n")
#     cat("\n")
#
#     contribs[,1] <- gdm.varcont(cdata, geo = geo)
#
#     # calculating variable contributions of permuted models
#
#     cat("Calculating variable contributions of permuted models\n")
#     cat("\n")
#
#     for (t in 1:w) {
#       cat(paste("Permutation",t,"of",w,"\n", sep = " "))
#       cat("\n")
#
#       q <- t+1
#       dist <- vegdist(as.data.frame(bio.perm$perm[t]), method=metric)
#       pdata <- cdata
#       pdata[,1] <- dist
#
#       contribs[,q] <- gdm.varcont(pdata, geo = geo)
#     }
#
#     # Feeding variable contribution matrix
#
#     cont.perm <- matrix(0,l1,w)
#
#     for (g in 1:l1){
#       for (f0 in 1:w){
#         f <- f0+1
#         cont.perm[g,f0] <- contribs[g,1] < contribs[g,f]
#       }
#     }
#
#     cont.sig <- matrix(0,l1,1)
#
#     # significance test
#
#     for (e in 1:l1){
#       cont.sig[e,] <- sum(cont.perm[e,])/w
#     }
#
#     cat("Significance test of GDM model variable contributions completed\n")
#     cat("\n")
#
#     return(cont.sig<=sig)
#   }
