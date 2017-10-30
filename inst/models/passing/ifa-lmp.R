################################################################################
## Test estimation of LMP model from rpf package
## Author: Carl F. Falk

require(OpenMx)
require(rpf)

data(LSAT6)
dat<-expandDataFrame(LSAT6, freqName="Freq")
dat<-mxFactor(as.data.frame(dat), levels=0:1)
N<-nrow(dat)
ni<-ncol(dat)
colnames(dat)<-paste("X",1:ni,sep="")


################################################################################
## k=0 is 1st order polynomial. Should be the same fit as 2PL

k<-0
par<-c(.1,0)
spec<-list()
spec[1:ni]<- list(rpf.lmp(k=k))

startingValues<-matrix(par,ncol=length(spec),
  nrow=length(par))

imat<-mxMatrix(name="item", values=startingValues, free=TRUE,
  dimnames=list(c("omega","xi"),colnames(dat)))

k0Model <- mxModel(model="k0Model", imat,
  mxData(observed=dat, type="raw"),
  mxExpectationBA81(spec,qwidth=5,qpoints=49),
    mxFitFunctionML(),
    mxComputeSequence(list(
      mxComputeEM('expectation', 'scores',
        mxComputeNewtonRaphson())
  )))

k0Model<-mxRun(k0Model)


## Try 2PL Model
par<-c(1,0)
spec<-list()
spec[1:ni]<-list(rpf.grm(2))

startingValues<-matrix(par,ncol=length(spec),
  nrow=length(par))

imat<-mxMatrix(name="item", values=startingValues, free=TRUE,
  dimnames=list(c("Slope","Intercept"),colnames(dat)))

twoplModel <- mxModel(model="twoplModel", imat,
  mxData(observed=dat, type="raw"),
  mxExpectationBA81(spec,qwidth=5,qpoints=49),
  mxFitFunctionML(),
  mxComputeSequence(list(
    mxComputeEM('expectation', 'scores',
      mxComputeNewtonRaphson())
  )))

twoplModel<-mxRun(twoplModel)

## 4933.307  
omxCheckCloseEnough(k0Model$output$fit,twoplModel$output$fit,.0001)

################################################################################
## k=1 is 3rd order polynomial. Test w/ priors vs. hardcoded param estimates
## obtained from CFF's own estimation code

k<-1
par<-c(.1,0,0,-1)

## prior means and variance
pmalpha<-0
pmtau<- -1
pvar<-50

spec<-list()
spec[1:ni]<-list(rpf.lmp(k=k))

startingValues<-matrix(par,ncol=length(spec),
  nrow=length(par))

imat<-mxMatrix(name="item", values=startingValues, free=TRUE,
  dimnames=list(c("omega","xi","alpha","tau"),colnames(dat)))

imat$labels['alpha',]<-paste('alpha',1:ni,sep="")
imat$labels['tau',]<-paste('tau',1:ni,sep="")

## Set up priors for alpha and tau; these are intended to be diffuse priors for stabilizing estmiation
## e.g., pi(tau) ~ N(-1,50), pi(alpha) ~ N(0,50), or with "50" replaced with some large value
## TO DO: use ifatools for gaussian priors
gaussPriorAlpha <- mxMatrix(name="gaussPriorAlpha",nrow=1, ncol=length(1:ni),
  free=TRUE, labels=imat$labels['alpha',1:ni],values=imat$values['alpha',1:ni])
gaussPriorTau <- mxMatrix(name="gaussPriorTau",nrow=1, ncol=length(1:ni),
  free=TRUE, labels=imat$labels['tau',1:ni],values=imat$values['tau',1:ni])

gaussMalpha <- mxMatrix(name="gaussMalpha", nrow=1, ncol=length(1:ni), values=pmalpha)
gaussSDalpha <- mxMatrix(name="gaussSDalpha", nrow=1, ncol=length(1:ni),values=sqrt(pvar))
gaussMtau <- mxMatrix(name="gaussMtau", nrow=1, ncol=length(1:ni), values=pmtau)
gaussSDtau <- mxMatrix(name="gaussSDtau", nrow=1, ncol=length(1:ni), values=sqrt(pvar))

## How much -2logLik changes due to the priors and resulting derivatives
gaussFitAlpha <- mxAlgebra(sum(log(2*pi) + 2*log(gaussSDalpha) +
  (gaussPriorAlpha-gaussMalpha)^2/gaussSDalpha^2), name="gaussFitAlpha")
gaussGradAlpha <- mxAlgebra(2*(gaussPriorAlpha - gaussMalpha)/gaussSDalpha^2, name="gaussGradAlpha",
  dimnames=list(c(),gaussPriorAlpha$labels))
gaussHessAlpha <- mxAlgebra(vec2diag(2/gaussSDalpha^2), name="gaussHessAlpha",
  dimnames=list(gaussPriorAlpha$labels, gaussPriorAlpha$labels))

gaussFitTau <- mxAlgebra(sum(log(2*pi) + 2*log(gaussSDtau) +
  (gaussPriorTau-gaussMtau)^2/gaussSDtau^2), name="gaussFitTau")
gaussGradTau <- mxAlgebra(2*(gaussPriorTau - gaussMtau)/gaussSDtau^2, name="gaussGradTau",
  dimnames=list(c(),gaussPriorTau$labels))
gaussHessTau <- mxAlgebra(vec2diag(2/gaussSDtau^2), name="gaussHessTau",
  dimnames=list(gaussPriorTau$labels, gaussPriorTau$labels))

gaussModelAlpha <- mxModel(model="gaussModelAlpha", gaussPriorAlpha, gaussMalpha, gaussSDalpha,
  gaussFitAlpha, gaussGradAlpha, gaussHessAlpha,
  mxFitFunctionAlgebra("gaussFitAlpha", gradient="gaussGradAlpha", hessian="gaussHessAlpha"))

gaussModelTau <- mxModel(model="gaussModelTau", gaussPriorTau, gaussMtau, gaussSDtau,
  gaussFitTau, gaussGradTau, gaussHessTau,
  mxFitFunctionAlgebra("gaussFitTau", gradient="gaussGradTau", hessian="gaussHessTau"))

itemModel <- mxModel(model="itemModel", imat,
  mxData(observed=dat, type="raw"),
  mxExpectationBA81(spec,qwidth=5,qpoints=49),
  mxFitFunctionML())

fit<-mxFitFunctionMultigroup(groups=c('itemModel.fitfunction',
  'gaussModelAlpha.fitfunction','gaussModelTau.fitfunction'))

k1PriorModel <- mxModel(model="k1PriorModel", itemModel, gaussModelAlpha, gaussModelTau,
  fit,
  mxComputeSequence(list(
    mxComputeEM('itemModel.expectation', 'scores',
      mxComputeNewtonRaphson(maxIter=10000),
      maxIter=2000)
                   )))

k1PriorModel<-mxRun(k1PriorModel)
  
## log likelihood if Bayesian estimates substituted into ML fit function
omxCheckCloseEnough(k1PriorModel$output$algebras$itemModel.fitfunction,4929.897,1e-3)

est<-k1PriorModel$output$estimate

myEst<-c(-0.27386336, 2.69911765,-0.06118075,-3.03394846,
  -0.10517433, 0.66718611,-0.72263680,-2.99031027,
  -0.60439972, 0.65458936, 1.23987805,-1.55414690,
  -2.09965081, 1.09230606,-2.87804350,-0.82664116,
  -0.91245939, 1.74431613,-1.42863857,-2.27822768)
  
## Should be close, but maybe not exact due to the scaling of some parameters
omxCheckCloseEnough(est,myEst,.1)
