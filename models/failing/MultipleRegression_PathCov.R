require(OpenMx)
myRegDataCov<-matrix(
	c(1.116, 0.539, 0.289, 
	  0.539, 0.933, 0.312,
	  0.289, 0.313, 0.836),
	nrow=3,
	dimnames=list(
	c("x","y","z"),c("x","y","z"))
	)

myRegDataMeans<-c(0.054, 0.574, 4.061)

model<-mxModel("Multiple Regression - Path", 
      type="RAM",
      mxData(myRegDataCov, type="cov", mean=myRegDataMeans, numObs=100),
      manifestVars=c("x","y","z"),
      mxPath(from="y", 
            arrows=2,
            free=TRUE, 
            values = 1,
            labels=c("residual")), # residual variances
      mxPath(from=c("x","z"), 
            arrows=2,
            free=TRUE, 
            all=TRUE,
            values = c(1,.5,.5,1),
            labels=c("varx","covxz","covxz","varz")), # exogenous covariance matrix
      mxPath(from=c("x","z"),
            to="y",
            arrows=1,
            free=TRUE,
            values=c(1,1),
            label=c("betax","betaz")), # regression weights
      mxPath(from="one",
            to=c("x","y","z"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1),
            labels=c("meanx","beta0","meanz")) # means
      ) # close model
      
multipleRegPathCov<-mxRun(model)

multipleRegPathCov@output

omxCheckCloseEnough(multipleRegPathCov@output$estimate[["beta0"]], 1.6331, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["betax"]], 0.4246, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["betaz"]], 0.2260, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["residual"]], 0.646, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["varx"]], 1.116, 0.001) multipleRegPathRaw omxCheckCloseEnough(multipleRegPathCov@output$estimate[["varz"]], 0.836, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["cov"]], 0.289, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["meanx"]], 0.054, 0.001)
omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["meanz"]], 4.061, 0.001)
