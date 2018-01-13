library(OpenMx)
library(MASS)

set.seed(1)

acePow2 <- function(add, com, Nmz, Ndz){
  
  nv <- 1
  ntv <- nv*2
  
  AA <- add
  CC <- com
  EE <- 1 - add - com
  
  mzMat <- matrix(c(AA + CC + EE, AA + CC, AA + CC, AA + CC + EE),2)
  dzMat <- matrix(c(AA + CC + EE, .5*AA + CC, .5*AA + CC, AA + CC + EE),2)
  
  mzData <- mvrnorm(Nmz , mu = c(0,0), mzMat, empirical = T)
  dzData <- mvrnorm(Ndz , mu = c(0,0), dzMat, empirical = T)
  
  selVars <- paste("t", 1:2, sep = "")
  colnames(mzData) <- colnames(dzData) <- selVars
  
  MZdata  <-  mxData(mzData, type="raw" )
  DZdata  <-  mxData(dzData, type="raw" )
  
  
  A <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=.6, label="A11", name="A" )
  C <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=.6, label="C11", name="C" )
  E <-  mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=.6, label="E11", name="E" )
  
  Mean    <-  mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= 0, label="mean", name="Mean" )
  expMean <-  mxAlgebra( expression= cbind(Mean,Mean), name="expMean")
  
  expCovMZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E , A+C),cbind(A+C   , A+C+E)), name="expCovMZ" )
  expCovDZ <- mxAlgebra( expression= rbind  ( cbind(A+C+E     , 0.5%x%A+C),cbind(0.5%x%A+C , A+C+E)),  name="expCovDZ" ) 
  
  obs <-  list(A,C,E,Mean,expMean,expCovMZ,expCovDZ)
  
  fun <- mxFitFunctionML()
  mzExp <- mxExpectationNormal(covariance="expCovMZ", means="expMean", dimnames=selVars )
  dzExp <- mxExpectationNormal(covariance="expCovDZ", means="expMean", dimnames=selVars )
  
  MZ <- mxModel("MZ", obs, MZdata, fun, mzExp)
  DZ <- mxModel("DZ", obs, DZdata, fun, dzExp)
  
  aceFun <- mxFitFunctionMultigroup(c("MZ","DZ"))
  ace <- mxModel("ACE", MZ, DZ, fun, aceFun)
  
  aceFit <- suppressWarnings(mxRun(ace, silent = T))
}

modA2 <- acePow2(add = .33, com = .3, Nmz = 1000, Ndz = 1000)

ceFit2 <- omxSetParameters(modA2, labels = "A11", free = F, values = 0)
ceFit2 <- mxRun(ceFit2)

got <- mxPowerSearch(modA2, ceFit2, probes = 50)
got <- mxPowerSearch(modA2, ceFit2, previousRun = got)

got2 <- mxPowerSearch(modA2, ceFit2, method = "ncp")

omxCheckCloseEnough(c(pmin(got2[,'power'] - got[,'lower'], 0),
                      pmin(got[,'upper'] - got2[,'power'], 0)),
                    rep(0,40), .01)

omxCheckEquals(mxPower(modA2, ceFit2, method='ncp'), 224)

omxCheckCloseEnough(mxPower(modA2, ceFit2, method='ncp', n=224, power=NULL), .817, .01)

omxCheckCloseEnough(mxPower(modA2, ceFit2, n=224, power=NULL, probes = 100),
                    .8, .5)

omxCheckCloseEnough(mxPower(modA2, ceFit2), 211.59, 5)
