# cdata.reduce <-
#   function(cdata,               # compiled dataset as output from "data.read" function
#            sigtest)             # predictor variable contribution significance test, as output by gdm.varsig
#   {
#
#     #
#     # p. j. leitao - 12th May 2014
#     #
#     # function for reducing compiled dataset based on significance test of GDM predictor variable contributions
#     #
#     # delivers reduced compiled dataset
#     #
#
#     # Reducing dataset
#
#     cat("Reducing compiled dataset following significance test result\n")
#     cat("\n")
#
#     cdata.sigtest <- rbind(as.matrix(as.logical(c("T","T","T","T","T","T"))),sigtest,sigtest)
#     new.cdata <- cdata[,cdata.sigtest]
#
#     cat("Data reduction of compiled dataset finished\n")
#     cat("\n")
#
#     return(new.cdata)
#   }
