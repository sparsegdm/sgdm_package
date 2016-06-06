"gdm.validate" <-
  function(caldata,                   # calibration data
           valdata,                   # validation data
           geo = F)                   # optional use of geographical distance as predictor in GDM model set as FALSE per default
  {

    # fitting  model

    model <- gdm(caldata, geo = geo)

    # calculating predicted dissimilarities for partial model

    predicted <- predict.gdm(model,valdata)
    predicted <- as.data.frame(predicted)

    # calculating model performance

    cat("\n")
    cat("Calculating model R-square (R2)...\n")
    cat("\n")
    performance <- (cor(valdata[,1], predicted)^2)*100

    cat("Model performance calculated\n")
    cat("\n")

    return(performance)
  }

