myRegDataCov<-matrix(
	c(1.116, 0.539,
	  0.539, 0.933),
	nrow=2,
	dimnames=list(
	c("x","y"),c("x","y"))
	)
	
myRegDataMeans<-c(0.05416, 2.57393)

model<-mxModel("Simple Regression, Path Specification", 
      type="RAM",
      mxData(myRegDataCov, type="cov", means=myRegDataMeans, numObs=100),
      manifestVars=c("x","y"),
      mxPath(from=c("x","y"), 
            arrows=2,
            free=TRUE, 
            values = c(1,1),
            labels=c("varx","residual")), # variances
      mxPath(from="x",
            to="y",
            arrows=1,
            free=TRUE,
            values=1,
            label="beta1"), # regression weight
      mxPath(from="one",
            to=c("x","y"),
            arrows=1,
            free=TRUE,
            values=c(1,1),
            labels=c("meanx","beta0")) # means
      ) # close model
      
regressionPathCov<-mxRun(model)

regressionPathCov@output