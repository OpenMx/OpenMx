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
	
myRegDataMeans<-c(2.582, 0.054, 2.574. 4.061)
	
model<-mxModel("Multivariate Regression (Path Analysis)- Path", 
      type="RAM",
      mxData(myRegDataCov,type="cov", mean=myRegDataMeans, numObs=100),
      manifestVars=c("w","x","y","z"),
      mxPath(from=c("w","y"), 
            arrows=2,
            free=TRUE, 
            all=FALSE,
            values = c(1,1),
            labels=c("residualw","residualy")), # residual variances
      mxPath(from=c("x","z"), 
            arrows=2,
            free=TRUE, 
            all=TRUE,
            values = c(1,.5,.5,1),
            labels=c("varx","cov","cov","varz")), # exogenous covariance matrix
      mxPath(from=c("x","z"),
            to="y",
            arrows=1,
            free=TRUE,
            values=c(1,1),
            label=c("betayx","betayz")), # regression weights y
      mxPath(from=c("x","z"),
            to="w",
            arrows=1,
            free=TRUE,
            values=c(1,1),
            label=c("betawx","betawz")), # regression weights w
      mxPath(from="one",
            to=c("w","x","y","z"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1,1),
            labels=c("betaw","meanx","betay","meanz")) # means
      ) # close model
      
multiRegPathCov<-mxRun(model)

multiRegPathCov@output

omxCheckCloseEnough(multiRegPathCov@output$estimate[["betay"]], 1.6331, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["betayx"]], 0.4246, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["betayz"]], 0.2260, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["residualy"]], 0.646, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["betaw"]], 0.51391, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["betawx"]], -0.23102, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["betawz"]], 0.51223, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["residualw"]], 0.60964, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["varx"]], 1.116, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["varz"]], 0.836, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["covxz"]], 0.289, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["meanx"]], 0.054, 0.001)
omxCheckCloseEnough(multiRegPathCov@output$estimate[["meanz"]], 4.061, 0.001)
