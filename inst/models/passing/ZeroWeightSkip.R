library(OpenMx)

set.seed(1)
data(myFADataRaw)

myFADataRaw <- myFADataRaw[,c("x1","x2","x3","x4","x5","x6")]
myFADataRaw$weight <- as.numeric(sample(0:2, nrow(myFADataRaw), replace=TRUE))

dataRaw      <- mxData( observed=myFADataRaw, type="raw", weight="weight" )
resVars      <- mxPath( from=c("x1","x2","x3","x4","x5","x6"), arrows=2,
                        free=TRUE, values=c(1,1,1,1,1,1),
                        labels=c("e1","e2","e3","e4","e5","e6") ) 
latVar       <- mxPath( from="F1", arrows=2,
                        free=TRUE, values=1, labels ="varF1" )
facLoads     <- mxPath( from="F1", to=c("x1","x2","x3","x4","x5","x6"), arrows=1,
                        free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=c(1,1,1,1,1,1),
                        labels =c("l1","l2","l3","l4","l5","l6") )
means        <- mxPath( from="one", to=c("x1","x2","x3","x4","x5","x6","F1"), arrows=1,
                        free=c(F,T,T,T,T,T,FALSE), values=c(0,1,1,1,1,1,0),
                        labels =c("rowMean1[1,1]","meanx2","meanx3",
                                  "meanx4","meanx5","meanx6",NA) ) 
means2       <- mxPath( from="one", to=c("x1","x2","x3","x4","x5","x6","F1"), arrows=1,
                        free=c(F,T,T,T,T,T,FALSE), values=c(0,1,1,1,1,1,0),
                        labels =c("meanx1","meanx2","meanx3",
                                  "meanx4","meanx5","meanx6",NA) ) 


oneFactorModel <- mxModel(
  "zws", type="RAM",
  manifestVars=c("x1","x2","x3","x4","x5","x6"), latentVars="F1",
  dataRaw, resVars, latVar, facLoads, means,
  mxMatrix("Full", 3, 1, values=c(-Inf,0,0), name="dmap"),
  mxAlgebra(dmap[1+data.weight,1], name="rowMean1"),
  mxComputeOnce('fitfunction', 'fit'))

oneFactorModel$fitfunction$jointConditionOn <- 'continuous'
f1 <- mxRun(oneFactorModel)      

oneFactorModel$fitfunction$jointConditionOn <- 'ordinal'
f2 <- mxRun(oneFactorModel)      

omxCheckCloseEnough(f1$output$fit, f2$output$fit, 1e-6)

omxCheckCloseEnough(f1$output$fit, 10013.433, .01)


#------------------------------------------------------------------------------
# Check non-integer weights work

oneFactorModel <- mxModel(
  "zws", type="RAM",
  manifestVars=c("x1","x2","x3","x4","x5","x6"), latentVars="F1",
  dataRaw, resVars, latVar, facLoads, means2,
  mxComputeOnce('fitfunction', 'fit'))

myFADataRaw$weight <- rep(0.5, nrow(myFADataRaw))
dataRaw      <- mxData( observed=myFADataRaw, type="raw", weight="weight" )
oneFactorModel <- mxModel(oneFactorModel, dataRaw)

f3 <- mxRun(oneFactorModel)

dataRaw      <- mxData( observed=myFADataRaw, type="raw")
oneFactorModel <- mxModel(oneFactorModel, dataRaw)

f4 <- mxRun(oneFactorModel)

omxCheckCloseEnough(logLik(f3), logLik(f4)/2, 1e-9)


