require(OpenMx)
myDataRaw<-read.table("myAutoregressiveData.txt",header=T)

myDataCov<-matrix(
	c(0.672, 0.315, 0.097, -0.037, 0.046,
	  0.315, 1.300, 0.428,  0.227, 0.146,
	  0.097, 0.428, 1.177,  0.568, 0.429,
	 -0.037, 0.227, 0.568,  1.069, 0.468,
	  0.046, 0.146, 0.429,  0.468, 1.031),
	nrow=5,
	dimnames=list(
	c("x1","x2","x3","x4","x5"),
	c("x1","x2","x3","x4","x5"))
	)
	
myDataMeans<-c(3.054, 1.385, 0.680, 0.254, -0.027)

model<-mxModel("Autoregressive Model - Path", 
      type="RAM",
      mxData(myDataCov,type="cov", means=myDataMeans, numObs=100),
      manifestVars=c("x1","x2","x3","x4","x5"),
      mxPath(from=c("x2","x3","x4","x5"),
            to=c("x1","x2","x3","x4"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1,1),
            labels=c("beta","beta","beta","beta")
            ),
      mxPath(from=c("x1","x2","x3","x4","x5"),
            arrows=2,
            free=TRUE,
            values=c(1,1,1,1,1),
            labels=c("varx1","e2","e3","e4","e5")
            ),
      mxPath(from="one",
            to=c("x1","x2","x3","x4","x5"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1,1,1),
            labels=c("mean1","mean2","mean.","mean4","mean5")
            )
      ) # close model
       
autoregressivePathCov<-mxRun(model)

autoregressivePathCov@output
