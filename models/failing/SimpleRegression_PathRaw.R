myRegDataRaw<-read.table("myRegData.txt",header=TRUE)

myRegDataRaw<-myRegDataRaw[,c("x","y")]

model<-mxModel("Simple Regression, Path Specification", 
      type="RAM",
      mxData(myRegDataRaw, type="raw"),
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
      
regressionPathRaw<-mxRun(model)

regressionPathRaw@output

#In model 'Simple Regression, Path Specification' NPSOL returned a non-zero status code 1. The final iterate x satisfies the optimality conditions to the accuracy requested, but the sequence of iterates has not yet converged. NPSOL was terminated because no further improvement could be made in the merit function. 

omxCheckCloseEnough(regressionPathRaw@output$estimate[["beta0"]], 2.54776, 0.001)
omxCheckCloseEnough(regressionPathRaw@output$estimate[["beta1"]], 0.48312, 0.001)
omxCheckCloseEnough(regressionPathRaw@output$estimate[["residual"]], 0.672, 0.01)
omxCheckCloseEnough(regressionPathRaw@output$estimate[["meanx"]], 0.05412, 0.001)
omxCheckCloseEnough(regressionPathRaw@output$estimate[["varx"]], 1.11654, 0.001)
