require(OpenMx)
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

oneFactorModel<-mxModel("Common Factor Model - Path", 
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
      

oneFactorFit<-mxRun(oneFactorModel)

oneFactorFit@output$estimate

omxCheckCloseEnough(oneFactorFit@output$estimate[["l2"]], 0.999, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l3"]], 0.959, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l4"]], 1.028, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l5"]], 1.008, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l6"]], 1.021, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["varF1"]], 0.645, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e1"]], 0.350, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e2"]], 0.379, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e3"]], 0.389, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e4"]], 0.320, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e5"]], 0.370, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e6"]], 0.346, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx1"]], 2.988, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx2"]], 3.011, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx3"]], 2.986, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx4"]], 3.053, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx5"]], 3.016, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx6"]], 3.010, 0.01)
