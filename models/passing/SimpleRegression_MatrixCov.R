require(OpenMx)
myRegDataCov<-matrix(
	c(1.116, 0.539,
	  0.539, 0.933),
	nrow=2,
	dimnames=list(
	c("x","y"),c("x","y"))
	)
	
myRegDataMeans<-c(0.05416, 2.57393)

model<-mxModel("Simple Regression - Matrix Specification", 
      mxData(myRegDataCov, type="cov", means=myRegDataMeans, numObs=100),
      mxMatrix("Full", nrow=2, ncol=2,
            values=c(0, 0,
                     1, 0),
            free=c(FALSE, FALSE,
                   TRUE,  FALSE),
            labels=c(NA,     NA,
                    "beta1", NA),
            byrow=TRUE,
            name="A"),
      mxMatrix("Symm", nrow=2, ncol=2, 
            values=c(1, 0,
                     0, 1),
            free=c(TRUE,  FALSE,
                   FALSE, TRUE),
            labels=c("varx", NA,
                      NA,    "residual"),
            byrow=TRUE,
            name="S"),
      mxMatrix("Iden",  nrow=2, ncol=2,
            name="F"),
      mxMatrix("Full", nrow=1, ncol=2,
            values=c(0,0),
            free=c(TRUE,TRUE),
            labels=c("meanx","beta0"),
            name="M"),
      mxRAMObjective("A","S","F","M")
      )
      
      
regressionMatrixCov<-mxRun(model)

regressionMatrixCov@output


#omxCheckCloseEnough(regressionMatrixCov@output$estimate[["beta0"]], 2.54776, 0.001)
omxCheckCloseEnough(regressionMatrixCov@output$estimate[["beta1"]], 0.48312, 0.001)
omxCheckCloseEnough(regressionMatrixCov@output$estimate[["residual"]], 0.672, 0.01)
#omxCheckCloseEnough(regressionMatrixRaw@output$estimate[["meanx"]], 0.05412, 0.001)
omxCheckCloseEnough(regressionMatrixCov@output$estimate[["varx"]], 1.11654, 0.001)