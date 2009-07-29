require(OpenMx)
myDataRaw<-read.table("myAutoregressiveData.txt",header=T)

model<-mxModel("Autoregressive Model, Matrix Specification, Raw Data",
      mxData(myDataRaw,type="raw"),
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
      
autoregressiveMatrixRaw<-mxRun(model)

#Comparing to Mplus values
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["beta"]], 0.427, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["varx"]], 0.665, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["e2"]], 1.142, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["e3"]], 1.038, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["e4"]], 0.791, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["e5"]], 0.818, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["mean1"]], 3.054, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["int2"]], 0.082, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["int3"]], 0.089, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["int4"]], -0.036, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw@output$estimate[["int5"]], -0.135, 0.001)
