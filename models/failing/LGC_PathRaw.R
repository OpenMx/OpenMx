require(OpenMx)
myLongitudinalData<-read.table("myLongitudinalData.txt",header=T)

model<-mxModel("Linear Growth Curve Model, Path Specification", 
      type="RAM",
      mxData(myLongitudinalData, type="raw"),
      manifestVars=c("x1","x2","x3","x4","x5"),
      latentVars=c("intercept","slope"),
      mxPath(from=c("x1","x2","x3","x4","x5"), 
            arrows=2,
            free=TRUE, 
            labels=c("residual","residual","residual","residual","residual"),
            values = c(1,1,1,1,1)), # constant residual variances
      mxPath(from=c("intercept","slope"), 
            arrows=2,
            free=TRUE, 
            labels=c("vari","vars"),
            values=c(1,1)), # latent variances
      mxPath(from="intercept", 
            to="slope", 
            arrows=2,
            free=TRUE, 
            labels="covis",
            values=.5), # latent covariance
      mxPath(from="intercept",
            to=c("x1","x2","x3","x4","x5"),
            arrows=1,
            free=FALSE,
            values=c(1,1,1,1,1)), # intercept loadings
      mxPath(from="slope",
            to=c("x1","x2","x3","x4","x5"),
            arrows=1,
            free=FALSE,
            values=c(0,1,2,3,4)), # slope loadings
      mxPath(from="one",
            to=c("x1","x2","x3","x4","x5"),
            arrows=1,
            free=FALSE,
            values=c(0,0,0,0,0)), # manifest means
      mxPath(from="one",
            to=c("intercept","slope"),
            arrows=1,
            free=TRUE,
            values=c(1,1),
            labels=c("meani","means")) # latent means
      ) # close model
      
growthPathRaw<-mxRun(model)
