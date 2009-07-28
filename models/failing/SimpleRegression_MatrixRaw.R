require(OpenMx)
myRegDataRaw<-read.table("myRegData.txt",header=TRUE)

myRegDataRaw<-myRegDataRaw[,c("x","y")]

model<-mxModel("Simple Regression - Matrix Specification", 
      mxData(myRegDataRaw,type="raw"),
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
      		dimnames=list(
      			c("x","y"),
      			c("x","y")),
            name="F"),
      mxMatrix("Full", nrow=1, ncol=2,
            values=c(0,0),
            free=c(TRUE,TRUE),
            labels=c("meanx","beta0"),
            dimnames=list(
      			NULL,
      			c("x","y")),
            name="M"),
      mxRAMObjective("A","S","F","M")
      )
      
regressionMatrixRaw<-mxRun(model)

regressionMatrixRaw@output


omxCheckCloseEnough(regressionMatrixRaw@output$estimate[["beta0"]], 2.54776, 0.001)
omxCheckCloseEnough(regressionMatrixRaw@output$estimate[["beta1"]], 0.48312, 0.001)
omxCheckCloseEnough(regressionMatrixRaw@output$estimate[["residual"]], 0.672, 0.001)
omxCheckCloseEnough(regressionMatrixRaw@output$estimate[["meanx"]], 0.05412, 0.001)
omxCheckCloseEnough(regressionMatrixRaw@output$estimate[["varx"]], 1.11654, 0.001)
# All pass except varx: varx = 1.1053196