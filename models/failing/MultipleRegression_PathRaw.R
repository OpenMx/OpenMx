require(OpenMx)
myRegDataRaw<-read.table("myRegData.txt",header=TRUE)

myRegDataRaw<-myRegDataRaw[,c("x","y","z")]

model<-mxModel("Multiple Regression - Path", 
      type="RAM",
      mxData(myRegDataRaw,type="raw"),
      manifestVars=c("x","y","z"),
      mxPath(from="y", 
            arrows=2,
            free=TRUE, 
            values = 1,
            labels=c("residual")), # residual variances
      mxPath(from=c("x","z"), 
            arrows=2,
            free=TRUE, 
            all=TRUE,
            values = c(1,.5,.5,1),
            labels=c("varx","covxz","covxz","varz")), # exogenous covariance matrix
      mxPath(from=c("x","z"),
            to="y",
            arrows=1,
            free=TRUE,
            values=c(1,1),
            label=c("betax","betaz")), # regression weights
      mxPath(from="one",
            to=c("x","y","z"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1),
            labels=c("meanx","beta0","meanz")) # means
      ) # close model
      
multipleRegPathRaw<-mxRun(model)

multipleRegPathRaw@output

omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["beta0"]], 1.6331, 0.001)
omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["betax"]], 0.4246, 0.001)
omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["betaz"]], 0.2260, 0.001)
omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["residual"]], 0.646, 0.001)
omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["varx"]], 1.116, 0.001) 
omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["varz"]], 0.836, 0.001)
omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["cov"]], 0.289, 0.001)
omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["meanx"]], 0.054, 0.001)
omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["meanz"]], 4.061, 0.001)
