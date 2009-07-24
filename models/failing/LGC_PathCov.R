require(OpenMx)
myLongitudinalDataCov<-matrix(
	c(6.362, 4.344, 4.915,  5.045,  5.966,
	  4.344, 7.241, 5.825,  6.181,  7.252,
	  4.915, 5.825, 9.348,  7.727,  8.968,
	  5.045, 6.181, 7.727, 10.821, 10.135,
	  5.966, 7.252, 8.968, 10.135, 14.220),
	nrow=5,
	dimnames=list(
		c("x1","x2","x3","x4","x5"),c("x1","x2","x3","x4","x5"))
	)

myLongitudinalDataMean<-c(9.864, 11.812, 13.612, 15.317, 17.178)

model<-mxModel("Linear Growth Curve Model, Path Specification", 
      type="RAM",
      mxData(myLongitudinalDataCov, type="cov", mean=myLongitudinalDataMean, numObs=500),
      manifestVars=c("x1","x2","x3","x4","x5"),
      latentVars=c("intercept","slope"),
      mxPath(from=c("x1","x2","x3","x4","x5"), 
            arrows=2,
            free=TRUE, 
            labels=c("residual","residual","residual","residual","residual"),
            values = c(1,1,1,1,1)), # constant residual variances
      mxPath(from=c("intercept","slope"), 
            arrows=2,
            free=TRUE, 
            labels=c("vari","vars"),
            values=c(1,1)), # latent variances
      mxPath(from="intercept", 
            to="slope", 
            arrows=2,
            free=TRUE, 
            labels="covis",
            values=.5), # latent covariance
      mxPath(from="intercept",
            to=c("x1","x2","x3","x4","x5"),
            arrows=1,
            free=FALSE,
            values=c(1,1,1,1,1)), # intercept loadings
      mxPath(from="slope",
            to=c("x1","x2","x3","x4","x5"),
            arrows=1,
            free=FALSE,
            values=c(0,1,2,3,4)), # slope loadings
      mxPath(from="one",
            to=c("x1","x2","x3","x4","x5"),
            arrows=1,
            free=FALSE,
            values=c(0,0,0,0,0)), # manifest means
      mxPath(from="one",
            to=c("intercept","slope"),
            arrows=1,
            free=TRUE,
            values=c(1,1),
            labels=c("meani","means")) # latent means
      ) # close model
      
growthPathCov<-mxRun(model)
