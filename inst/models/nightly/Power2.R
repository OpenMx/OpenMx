library(OpenMx)
library(MASS)

suppressWarnings(RNGversion("3.5"))
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

# Build a true and a false model
ACEfit <- acePow2(add = .33, com = .3, Nmz = 1000, Ndz = 1000)
# drop A
CEfit <- omxSetParameters(ACEfit, labels = "A11", free = FALSE, values = 0)
CEfit <- mxRun(CEfit)

AEfit <- omxSetParameters(ACEfit, labels = "C11", free = FALSE, values = 0)
AEfit <- mxRun(AEfit)

# ACE and CE coefficients
round(coef(ACEfit),2)
round( coef(CEfit),2)
round( coef(AEfit),2)

# Check N estimated for 80% power is 0.223 pairs total
omxCheckEquals(as.numeric(mxPower(ACEfit, CEfit, method='ncp')), 446)

# Set N=500, request power and check it's within .01 of .8
omxCheckCloseEnough(as.numeric(mxPower(ACEfit, CEfit, method='ncp', n=500, power=NULL)), .857, .01)

# Set N=500, request power with method empirical, 100 probes, and check it's +/- 5% of .8
omxCheckCloseEnough(as.numeric(mxPower(ACEfit, CEfit, n=500, power=NULL, probes = 100)), .8, .5)

# Empirically search for N required to reject false model (CEfit) 80% of long-run occasions
got <- mxPower(ACEfit, CEfit)
omxCheckCloseEnough(as.numeric(got), 458)


got <- mxPowerSearch(ACEfit, CEfit, probes = 50)
got <- mxPowerSearch(ACEfit, CEfit, previousRun = got,
                     grid=seq(100,800,length.out = 20))

got2 <- mxPowerSearch(ACEfit, CEfit, method = "ncp",
                      grid=seq(100,800,length.out = 20))

omxCheckCloseEnough(c(pmin(got2[,'power'] - got[,'lower'], 0),
                      pmin(got[,'upper'] - got2[,'power'], 0)),
                    rep(0,40), .01)

