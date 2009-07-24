require(OpenMx)
myDataRaw<-read.table("myAutoregressiveData.txt",header=T)

model<-mxModel("Autoregressive Model - Path", 
      type="RAM",
      mxData(myDataRaw,type="raw"),
      #mxData(myDataCov,type="cov", means=myDataMeans, numObs=100),
      manifestVars=c("x1","x2","x3","x4","x5"),
      mxPath(from=c("x2","x3","x4","x5"),
            to=c("x1","x2","x3","x4"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1,1),
            labels=c("beta","beta","beta","beta")
            ),
      mxPath(from=c("x1","x2","x3","x4","x5"),
            arrows=2,
            free=TRUE,
            values=c(1,1,1,1,1),
            labels=c("varx1","e2","e3","e4","e5")
            ),
      mxPath(from="one",
            to=c("x1","x2","x3","x4","x5"),
            arrows=1,
            free=TRUE,
            values=c(1,1,1,1,1),
            labels=c("mean1","mean2","mean.","mean4","mean5")
            )
      ) # close model
       
autoregressivePathRaw <-mxRun(model)

autoregressivePathRaw@output
