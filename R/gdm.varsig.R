#' @title Calculates statistical significance of predictor variable contributions in GDM model
#'
#' @description
#' This function calculates the statistical significance of predictor variable contributions in GDM model, as calculated with the \code{gdm.varcont} function.
#'
#' Significance is calculated through random matrix permutations, to infer if the variable (drop) contributions are greater than those obtained by chance.
#'
#' It requires a biological ("bioData" format) and a predictor ("spData" format) datasets, the dissimilarity metric to be used, the optional use of geographical distance as predictor variable in the GDM, the number of permutations to be run and the significance level to be tested.
#'
#' This current implementation only allows biological data in the format 1 using abundance values, as described in the \code{gdm} package.
#'
#' For more details relating to "bioData" and "predData" data formats, check \code{gdm} package.
#'
#' @param predData Predictor dataset ("predData" format).
#' @param bioData Biological dataset ("bioData" format).
#' @param metric Dissimilarity metric to be used.
#' @param geo Optional use of geographical distance as predictor in GDM model. Set to \code{FALSE} per default
#' @param perm Number of matrix permutations to be used for calculating statistical significance.
#' @param sig Significance level (p-value) to be tested. Set as 0.05 per default.
#' @return Returns vector resulting of the significance test. Variables with significant model contributions are assigned value \code{TRUE}, and non-significant variables the value \code{FALSE}, according to the defined p-value threshold.
#' @export

gdm.varsig <-
  function(predData,
           bioData,
           metric = "bray",
           geo = F,
           perm = 100,
           sig = 0.05)
  {

    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("raster package needed for this function to work. Please install it.",
           call. = FALSE)
    }

    # data reading

    spData <- gdm::formatsitepair(bioData, 1, dist = "bray", abundance = TRUE,
                             siteColumn = "Plot_ID", XColumn = "X",YColumn = "Y",
                             predData = predData)

    l1 <- (ncol(spData)-6)/2
    j2 <- ncol(bioData)

    w <- perm
    w1 <- w+1

    cat("Testing for significance of GDM model variable contributions\n")
    cat("\n")


    # creating variable contribution matrix

    contribs <- matrix(0,l1,w1)

    # initiating significance test: species matrix permutations

    cat("Checking statistical significance of GDM variable contributions:\n")
    cat("Performing species matrix permutations\n")
    cat("\n")

    bio.m <- as.matrix(bioData[,2:j2])
    bio.perm <- vegan::permatfull(bio.m, fixedmar = "row", times = w)

    perm.m <- matrix(0,w,1)

    # calculating variable contributions of original model

    cat("Running original GDM model and calculating variable contributions\n")
    cat("\n")

    contribs[,1] <- gdm.varcont(spData, geo = geo)

    # calculating variable contributions of permuted models

    cat("Calculating variable contributions of permuted models\n")
    cat("\n")

    for (t in 1:w) {
      cat(paste("Permutation",t,"of",w,"\n", sep = " "))
      cat("\n")

      q <- t+1
      dist <- vegan::vegdist(as.data.frame(bio.perm$perm[t]), method=metric)
      pdata <- spData
      pdata[,1] <- dist

      contribs[,q] <- gdm.varcont(pdata, geo = geo)
    }

    # Feeding variable contribution matrix

    cont.perm <- matrix(0,l1,w)

    for (g in 1:l1){
      for (f0 in 1:w){
        f <- f0+1
        cont.perm[g,f0] <- contribs[g,1] < contribs[g,f]
      }
    }

    cont.sig <- matrix(0,l1,1)

    # significance test

    for (e in 1:l1){
      cont.sig[e,] <- sum(cont.perm[e,])/w
    }

    cat("Significance test of GDM model variable contributions completed\n")
    cat("\n")

    return(cont.sig<=sig)
  }
