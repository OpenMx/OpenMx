require(OpenMx)
myFADataCov<-matrix(
	c(0.997, 0.642, 0.611, 0.342, 0.299, 0.337,
	  0.642, 1.025, 0.608, 0.273, 0.282, 0.287,
	  0.611, 0.608, 0.984, 0.286, 0.287, 0.264,
	  0.342, 0.273, 0.286, 0.993, 0.472, 0.467,
	  0.299, 0.282, 0.287, 0.472, 0.978, 0.507,
	  0.337, 0.287, 0.264, 0.467, 0.507, 1.059),
	nrow=6,
	dimnames=list(
		c("x1","x2","x3","y1","y2","y3"),
		c("x1","x2","x3","y1","y2","y3"))
	)

myFADataMeans<-c(2.988, 3.011, 2.986, 2.955, 2.956, 2.967)

model<-mxModel("Two Factor Model - Path", 
      type="RAM",
      mxData(myFADataCov, type="cov", mean=myFADataMeans, numObs=500),
      manifestVars=c("x1","x2","x3","x4","x5","x6"),
      latentVars=c("F1","F2"),
      mxPath(from=c("x1","x2","x3","x4","x5","x6"),
            arrows=2,
            free=TRUE,
            values=c(1,1,1,1,1,1),
            names=c("e1","e2","e3","e4","e5","e6")
            ), # residual variances
      mxPath(from=c("F1","F2"),
            arrows=2,
            all=2,
            free=TRUE,
            values=c(1, .5,
                    .5, 1),
            names=c("varF1","cov","cov","varF2")
            ), # latent variance
      mxPath(from="F1",
            to=c("x1","x2","x3"),
            arrows=1,
            free=c(FALSE,TRUE,TRUE),
            values=c(1,1,1),
            names=c("l1","l2","l3")
            ), # factor loadings, factor 1
      mxPath(from="F2",
            to=c("x4","x5","x6"),
            arrows=1,
            free=c(FALSE,TRUE,TRUE),
            values=c(1,1,1),
            names=c("l4","l5","l6")
            ), # factor loadings, factor 1
      mxPath(from="one",
            to=c("x1","x2","x3","x4","x5","x6","F1","F2"),
            arrows=1,
            free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE),
            values=c(1,1,1,1,1,1,0,0),
            names=c("meanx1","meanx2","meanx3",
                    "meanx4","meanx5","meanx6",
                    NA,NA)
            ) # means
      ) # close model
      
 
twoFactorPathCov<-mxRun(model)
