myRegDataCov<-read.table("myRegData.txt",header=TRUE)
	
model<-mxModel("Multivariate Regression (Path Analysis)- Path", 
      type="RAM",
      mxData(myRegDataRaw,type="raw"),
      manifestVars=c("w","x","y","z"),
      mxPath(from=c("w","y"), 
            arrows=2,
            free=TRUE, 
            all=FALSE,
            values = c(1,1),
            labels=c("residualw","residualy")), # residual variances
      mxPath(from=c("x","z"), 
            arrows=2,
            free=TRUE, 
            all=TRUE,
            values = c(1,.5,.5,1),
            labels=c("varx","cov","cov","varz")), # exogenous covariance matrix
      mxPath(from=c("x","z"),
            to="y",
            arrows=1,
            free=TRUE,
            values=c(1,1),
            label=c("betayx","betayz")), # regression weights y
      mxPath(from=c("x","z"),
            to="w",
            arrows=1,
            free=TRUE,
            values=c(1,1),
            label=c("betawx","betawz")), # regression weights w
      mxPath(from="one",
            to=c("w","x","y","z"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1,1),
            labels=c("betaw","meanx","betay","meanz")) # means
      ) # close model
      
multiRegPathRaw<-mxRun(model)

multiRegPathRaw@output

omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betay"]], 1.6331, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betayx"]], 0.4246, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betayz"]], 0.2260, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["residualy"]], 0.646, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betaw"]], 0.51391, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betawx"]], -0.23102, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betawz"]], 0.51223, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["residualw"]], 0.60964, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["varx"]], 1.116, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["varz"]], 0.836, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["covxz"]], 0.289, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["meanx"]], 0.054, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["meanz"]], 4.061, 0.001)