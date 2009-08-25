require(OpenMx)
myRegDataCov<-matrix(
	c(0.808,-0.110, 0.089, 0.361,
	 -0.110, 1.116, 0.539, 0.289,
	  0.089, 0.539, 0.933, 0.312,
	  0.361, 0.289, 0.312, 0.836),
	nrow=4,
	dimnames=list(
		c("w","x","y","z"),
		c("w","x","y","z"))
	)
	
myRegDataMeans<-c(2.582, 0.054, 2.574, 4.061)
	
multivariateRegModel<-mxModel("Multiple Regression - Matrix Specification", 
      mxData(myRegDataCov,type="cov", mean=myRegDataMeans, numObs=100),
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
               F, F, T, F,
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
    mxMatrix("Full", nrow=4, ncol=1,
        values=c(0,0,0,0),
        free=c(T,T,T,T),
        labels=c("betaw","meanx","betay","meanz"),
        dimnames=list(
      	    c("w","x","y","z"), NULL),
        name="M"),
    mxRAMObjective("A","S","F","M")
)
      
multivariateRegFit<-mxRun(multivariateRegModel)

multivariateRegFit@output

# Old Mx Answers
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betay"]], 1.6312, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayx"]], 0.4243, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayz"]], 0.2265, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualy"]], 0.6336, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betaw"]], 0.5165, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawx"]], -0.2311, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawz"]], 0.5117, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualw"]], 0.5979, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["varx"]], 1.1160, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["varz"]], 0.8360, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["covxz"]], 0.2890, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanx"]], 0.0540, 0.001)
omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanz"]], 4.0610, 0.001)


# omxCheckCloseEnough(multivariateRegFit@output$estimate[["betay"]], 1.6331, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayx"]], 0.4246, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["betayz"]], 0.2260, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualy"]], 0.646, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["betaw"]], 0.51391, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawx"]], -0.23102, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["betawz"]], 0.51223, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["residualw"]], 0.60964, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["varx"]], 1.116, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["varz"]], 0.836, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["covxz"]], 0.289, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanx"]], 0.054, 0.001)
# omxCheckCloseEnough(multivariateRegFit@output$estimate[["meanz"]], 4.061, 0.001)
