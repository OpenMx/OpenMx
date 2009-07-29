require(OpenMx)
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

model<-mxModel("Autoregressive Model, Matrix Specification, Covariance Data",
      mxData(myDataCov,type="cov", means=myDataMeans, numObs=100),
      mxMatrix("Full", nrow=5, ncol=5,
            values=c(0,1,0,0,0,
                     0,0,1,0,0,
                     0,0,0,1,0,
                     0,0,0,0,1,
                     0,0,0,0,0),
            free=c(F, T, F, F, F,
                   F, F, T, F, F,
                   F, F, F, T, F,
                   F, F, F, F, T,
                   F, F, F, F, F),
            labels=c(NA, "beta", NA,    NA,    NA,
                     NA, NA,    "beta", NA,    NA,
                     NA, NA,     NA,   "beta", NA,
                     NA, NA,     NA,    NA,   "beta",
                     NA, NA,     NA,    NA,    NA),
            byrow=TRUE,
            name="A"),
      mxMatrix("Symm", nrow=5, ncol=5, 
            values=c(1, 0, 0, 0, 0,
                     0, 1, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 1, 0,
                     0, 0, 0, 0, 1),
            free=c(T, F, F, F, F,
                   F, T, F, F, F,
                   F, F, T, F, F,
                   F, F, F, T, F,
                   F, F, F, F, T),
            labels=c("varx", NA,  NA,  NA,  NA,
                     NA,     "e2", NA,  NA,  NA,
                     NA,      NA, "e3", NA,  NA,
                     NA,      NA,  NA, "e4", NA,
                     NA,      NA,  NA,  NA, "e5"),
            byrow=TRUE,
            name="S"),
      mxMatrix("Iden",  nrow=5, ncol=5,
      		dimnames=list(
				c("x1","x2","x3","x4","x5"), c("x1","x2","x3","x4","x5")),
            name="F"),
      mxMatrix("Full", nrow=1, ncol=5,
            values=c(1,1,1,1,1),
            free=c(T, T, T, T, T),
            labels=c("mean1","int2","int3","int4","int5"),
            dimnames=list(
				NULL, c("x1","x2","x3","x4","x5")),
            name="M"),
      mxRAMObjective("A","S","F","M")
      )
      
autoregressiveMatrixCov<-mxRun(model)

autoregressiveMatrixCov@output

#Comparing to Mplus values
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["beta"]], 0.427, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["varx"]], 0.665, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["e2"]], 1.142, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["e3"]], 1.038, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["e4"]], 0.791, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["e5"]], 0.818, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["mean1"]], 3.054, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["int2"]], 0.082, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["int3"]], 0.089, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["int4"]], -0.036, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["int5"]], -0.135, 0.001)

#Comparing to Raw Data values
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["beta"]], 0.37298, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["varx"]], 0.61242, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["e2"]], 1.13246, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["e3"]], 0.89346, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["e4"]], 0.85459, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["e5"]], 1.02045, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["mean1"]], 3.05369, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["int2"]], 1.38547, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["int3"]], 0.68021, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["int4"]], 0.25410, 0.001)
omxCheckCloseEnough(autoregressiveMatrixCov@output$estimate[["int5"]], -0.02715, 0.001)
