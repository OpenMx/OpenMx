myFADataCov<-matrix(
	c(0.997, 0.642, 0.611, 0.672, 0.637, 0.677,
	  0.642, 1.025, 0.608, 0.668, 0.643, 0.676,
	  0.611, 0.608, 0.984, 0.633, 0.657, 0.626,
	  0.672, 0.668, 0.633, 1.003, 0.676, 0.665,
	  0.637, 0.643, 0.657, 0.676, 1.028, 0.654,
	  0.677, 0.676, 0.626, 0.665, 0.654, 1.020),
	nrow=6,
	dimnames=list(
		c("x1","x2","x3","x4","x5","x6"),
		c("x1","x2","x3","x4","x5","x6"))
	)

myFADataMeans<-c(2.988, 3.011, 2.986, 3.053, 3.016, 3.010)

model<-mxModel("Common Factor Model - Path", 
      type="RAM",
      mxData(myFADataCov, type="cov", mean=myFADataMeans, numObs=500),
      manifestVars=c("x1","x2","x3","x4","x5","x6"),
      latentVars="F1",
      mxPath(from=c("x1","x2","x3","x4","x5","x6"),
            arrows=2,
            free=TRUE,
            values=c(1,1,1,1,1,1),
            labels=c("e1","e2","e3","e4","e5","e6")
            ), # residual variances
      mxPath(from="F1",
            arrows=2,
            free=TRUE,
            values=1,
            labels ="varF1"
            ), # latent variance
      mxPath(from="F1",
            to=c("x1","x2","x3","x4","x5","x6"),
            arrows=1,
            free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
            values=c(1,1,1,1,1,1),
            labels =c("l1","l2","l3","l4","l5","l6")
            ), # factor loadings
      mxPath(from="one",
            to=c("x1","x2","x3","x4","x5","x6","F1"),
            arrows=1,
            free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE),
            values=c(1,1,1,1,1,1,0),
            labels =c("meanx1","meanx2","meanx3",
                    "meanx4","meanx5","meanx6",
                    NA)
            ) # means
      ) # close model
      

oneFactorPathCov<-mxRun(model)

oneFactorPathCov@output$estimate
