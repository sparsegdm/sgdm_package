gdm.map <- function(spData,        # Site pair table as from Formatsitetable gdm function
                    predMap,       # Data frame that has the same predictior variables as used for gdm model building
                    model,         # gdm.model
                    output="m",    # type of output: "n" NMDS model object
                                   #                 "p" NMDS transformed prediction map as data frame
                                   #                 "m" NMDS transformed prediction map as raster object (provide R raster object)
                    rst,           # R raster object with extent of the prediction map
                    k=0,           # number of NMDS components to extract; if not specified number of components will be derived be NMDS stress value
                    t=0.1)         # NMDS stress value threshold to extract number of components if k is not specified
{
  require(vegan)
  require(gdm)
  require(raster)

  cat("\n")
  cat("Starting beta diversity map prediction.")
  cat("\n")

  # function to derive map a gdm model output as
  # a) a NMDS transformed data frame
  # b) a raster image


  # derive information about input data
  pairs <- nrow(spData)
  n2 <- (1+sqrt(1+8*pairs))/2
  nVars<-(ncol(spData)-6)/2
  dummy_ID<-data.frame(ID=1:n2)

  # derive envdata from ipnut table
  first_sample<-cbind(dummy_ID[1,1], spData[1,3:4], spData[1,7:(nVars+6)])
  envdata_ <- cbind(dummy_ID[2:n2,1],spData[1:(n2-1),5:6],spData[1:(n2-1),(nVars+7):ncol(spData)])
  names<-c("ID","X","Y", paste0("Pred.", 1:(nVars)))
  colnames(first_sample)<-names
  colnames(envdata_)<-names
  envdata<-rbind(first_sample, envdata_)
  #dummy <- matrix(data=0, nrow=nrow(envdata), ncol=1)


  # predict dissimilarities for sample pairs
  sample.pair <- spData
  sample.pair.diss<-predict.gdm(model, sample.pair)

  #sample.pair.diss.mat<-as.data.frame(sample.pair.diss)

  X <- diag(x=0, nrow(envdata),nrow(envdata))
  X[upper.tri(X, diag=FALSE)] <- sample.pair.diss
  X <- X + t(X) - diag(diag(X))

  sample.pair.diss.mat<-as.matrix(X)

  if (k == 0){
    cat("\n")
    cat("Derive number of NMDS components based on stress values below" ,t, "\n")
    cat("\n")

    # derive k for NMDS based on 20 iterations and all possible components
    stress<-matrix(nrow=20,ncol=nrow(sample.pair.diss.mat))
    stress<-as.data.frame(stress)

    for (i in 1:nrow(sample.pair.diss.mat)){
      for (j in 1:20){

        table_nmds <- monoMDS(sample.pair.diss.mat, k=i, model ="global", maxit=1000)

        # Stress values of model
        # stressplot(table_nmds)
        stress[j,i]<-table_nmds$stress

      }

      cat("Mean stress value after 20 iterations with k =",i, "is", apply(stress[i],2, mean),"\n")
      cat("\n")
    }

    mean_stress<-as.matrix(apply(stress,2, mean))
    plot(mean_stress, main="Mean NMDS stress values out of 20 iterations")
    k<- as.numeric(length(which(mean_stress>t)))


  }

  # Perform NMDS on sample dissimialirites

  cat("\n")
  cat("Performing NMDS transformation with",k, "components.")
  cat("\n")

  sample_nmds <- monoMDS(sample.pair.diss.mat, k=k, model="global", maxit = 1000)

  print(sample_nmds$stress)

  stressplot(sample_nmds)

  if(output=="n"){

    cat("\n")
    cat("Done with NMDS transformation of the predicted sample dissimilarities. Output type R monoMDS object.")
    cat("\n")

    return(sample_nmds)

  }


  # Extract scores from NMDS model
  nmds_scores <- as.matrix(scores(sample_nmds))

  ####### here could be a potential output if just nmds should be extracted

  # prepare sample/map site pair table

  # check if both input data fit each other
  cat("\n")
  cat("Checking consistency of environmental map and sample data: ")

  j1 <- ncol(envdata)
  j3 <- ncol(predMap)
  n.sites <- nrow(envdata)
  weight<-1

  if (j1 != j3){
    stop("Map and samples not consistent, check input data!\n")
  }

  cat("Data check OK!\n")
  cat("\n")

  cat("Compiling dataset...\n")
  cat("\n")

  # creat site pair table
  map <- vector(mode = "list", length = n.sites)
  map <- lapply(1:length(map), function(x)map[[x]] <- predMap)
  map <- do.call("rbind", map)

  sites <- vector(mode = "list", length = n.sites)
  sites <- lapply(1:n.sites, function(x)do.call("rbind", rep(list(envdata[x,]), nrow(predMap))))
  sites <- do.call("rbind", sites)

  map.sample.pair <- cbind(rep(1, nrow(map)), rep(weight, nrow(map)), map[,2:3], sites[,2:3], map[,4:ncol(map)], sites[,4:ncol(sites)])
  names(map.sample.pair) <- c("distance", "weight", "p.x", "p.y", "s.x", "s.y", paste0("p.env.", 1:(ncol(map)-3)), paste0("s.env.", 1:(ncol(sites)-3)))


  # predict dissimilarities for map sample site table based on gdm model
  map.sample.pair.diss <- predict.gdm(model, map.sample.pair)

  # prepare predictions for NMDS transformation
  prediction.map <- matrix(map.sample.pair.diss, ncol=nrow(envdata), byrow=TRUE)

  # transform map sample site table predictions to NMDS components
  map_trans <- prediction.map %*% nmds_scores

  if(output == "p"){

    cat("\n")
    cat("Done with dissimilaritiy prediction and NMDS transformation. Output type R dataFrame.")
    cat("\n")

    return(map_trans)

  }


  rst_1<-raster(rst, layer=1)
  rst_2<-raster(rst, layer=1)
  rst_3<-raster(rst, layer=1)
  nmds_1<-map_trans[,1]
  nmds_2<-map_trans[,2]
  nmds_3<-map_trans[,3]
  values(rst_1)<-nmds_1
  values(rst_2)<-nmds_2
  values(rst_3)<-nmds_3

  nmds_map<-stack(rst_1,rst_2,rst_3)

  plotRGB(nmds_map, r=1, g=2, b=3, stretch="hist")

  if(output == "m"){

    cat("\n")
    cat("Done with dissimilaritiy prediction and NMDS transformation. Output type R Raster object.")
    cat("\n")

    return(nmds_map)

  }

}
