require(OpenMx)
myRegDataRaw<-read.table("myRegData.txt",header=TRUE)
	
multivariateRegModel<-mxModel("Multiple Regression - Matrix Specification", 
    mxData(myRegDataRaw,type="raw"),
    mxMatrix("Full", nrow=4, ncol=4,
        values=c(0,1,0,1,
                 0,0,0,0,
                 0,1,0,1,
                 0,0,0,0),
        free=c(F, T, F, T,
               F, F, F, F,
               F, T, F, T,
               F, F, F, F),
        labels=c(NA, "betawx", NA, "betawz",
                 NA,  NA,     NA,  NA, 
                 NA, "betayx", NA, "betayz",
                 NA,  NA,     NA,  NA),
        byrow=TRUE,
        name="A"),
    mxMatrix("Symm", nrow=4, ncol=4, 
        values=c(1,  0, 0,  0,
                 0,  1, 0, .5,
                 0,  0, 1,  0,
                 0, .5, 0,  1),
        free=c(T, F, F, F,
               F, T, F, T,
               F, T, F, T),
        labels=c("residualw",  NA,     NA,         NA,
                  NA,         "varx",  NA,        "covxz",
                  NA,          NA,    "residualy", NA,
                  NA,         "covxz", NA,        "varz"),
        byrow=TRUE,
        name="S"),
    mxMatrix("Iden",  nrow=4, ncol=4,
        dimnames=list(
            c("w","x","y","z"),
            c("w","x","y","z")),
        name="F"),
    mxMatrix("Full", nrow=1, ncol=4,
        values=c(0,0,0,0),
        free=c(T,T,T,T),
        labels=c("betaw","meanx","betay","meanz"),
        dimnames=list(
      	    NULL,c("w","x","y","z")),
        name="M"),
    mxRAMObjective("A","S","F","M")
)
      
multivariateRegFit<-mxRun(multivariateRegModel)

multivariateRegFit@output

omxCheckCloseEnough(multivariateRegFit@output$estimate[["betay"]], 1.6331, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayx"]], 0.4246, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayz"]], 0.2260, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualy"]], 0.646, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betaw"]], 0.51391, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawx"]], -0.23102, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawz"]], 0.51223, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualw"]], 0.60964, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["varx"]], 1.116, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["varz"]], 0.836, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["covxz"]], 0.289, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanx"]], 0.054, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanz"]], 4.061, 0.001)
