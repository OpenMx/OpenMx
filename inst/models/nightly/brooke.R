# Test result of algebras that depend on definition variables
# after mxRun. Contributed by Gregory Carey.

# full text for "Sat.fitFuntion not finite" is
# Error: The job for model 'Sat' exited abnormally with the error message: MxComputeGradientDescent: fitfunction Sat.fitfunction is not finite ()
# In addition: Warning message:
#   In model 'Sat' Optimizer returned a non-zero status code 10. Starting values are not feasible. Consider mxTryHard() 
#set.seed(101276) # gives "Sat.fitFunction not finite" works w mod 1
#set.seed(50074)  # gives "Sat.fitFunction not finite" does not work w mod1
set.seed(120450) # gives "Sat.fitFunction not finite" This works w mod1
#set.seed(270156) # gives "Sat.fitFunction not finite"
h1sq <- .62
h2sq <- .47
rg <- .55
re <- .34
covG <- sqrt(h1sq)*sqrt(h2sq)*rg
rp <- covG + sqrt((1 - h1sq)*(1 - h2sq))*re
pmat <- matrix(c(1, rp, rp, 1), 2, 2)
gmat <- matrix(c(h1sq, covG, covG, h2sq), 2, 2)
SigmaMZ <- rbind(cbind(pmat, gmat),
                 cbind(gmat, pmat))
SigmaDZ <- rbind(cbind(pmat, .5*gmat),
                 cbind(.5*gmat, pmat))
thresh <- c(.3, .5, .3, .5)
N <- 250



fcn <- function(N, Sigma, thresh) {
  age <- round(25 + 5*rnorm(N))
  X <- mvtnorm::rmvnorm(N, rep(0, ncol(Sigma)), Sigma)
  for (i in 1:4) {
    one <- which(X[,i] > thresh[i])
    zero <- which(X[,i] <= thresh[i])
    X[one,i] <- 1
    X[zero,i] <- 0
  }
  return(cbind(X, age, age))
}

temp <- fcn(N, SigmaMZ, thresh)
df <- cbind(rep(1, N), temp)
temp <- fcn(N, SigmaMZ, thresh)
temp <- cbind(rep(2, N), temp)
df <- rbind(df, temp)
temp <- fcn(N, SigmaDZ, thresh)
temp <- cbind(rep(3, N), temp)
df <- rbind(df, temp)
temp <- fcn(N, SigmaDZ, thresh)
temp <- cbind(rep(4, N), temp)
df <- rbind(df, temp)
colnames(df) <- c("zyg", "mj1", "nssi1", "mj2",
                  "nssi2", "age1", "age2")
mj_nssi <- data.frame(df, stringsAsFactors=F)        

require(OpenMx)

# --- this from Brooke's post 
#     http://openmx.ssri.psu.edu/thread/3996
#     note that start values have been changed to
#     reflect the parameters for the simulated data

#Make Factors
mj_nssi$mj1 <-mxFactor(mj_nssi$mj1, levels=c(0:1), exclude = NA)
mj_nssi$mj2  <-mxFactor(mj_nssi$mj2, levels=c(0:1), exclude = NA)
mj_nssi$nssi1 <-mxFactor(mj_nssi$nssi1, levels=c(0:1), exclude = NA)
mj_nssi$nssi2  <-mxFactor(mj_nssi$nssi2, levels=c(0:1), exclude = NA)

#Selects Vars for Analysis
vars <-c('nssi','mj')
selVars	<-c('nssi1', 'mj1', 'nssi2', 'mj2' )
useVars	<-c('nssi1', 'mj1', 'nssi2', 'mj2', 'age1', 'age2')

# Select Data for Analysis
mzfData <- subset(mj_nssi, zyg==1, useVars)
mzmData <- subset(mj_nssi, zyg==2, useVars)
dzfData <- subset(mj_nssi, zyg==3, useVars)
dzmData <- subset(mj_nssi, zyg==4, useVars)

nv <- 2	# number of variables per twin
ntv <- nv*2	# number of variables per pair
nth	<- 1	# number of max thresholds

# 1) Fits a constrained Polychoric correlation model
# TH same across twins but different across zyg groups
# Age effect is different across variables, but same across thresholds within variables (if c>2) ###?
# There is one overall rPH between var1-2 and the x-trait x-twin correlations are symmetric
# ------------------------------------------------------------------------------------------------------------------------------

# CREATE LABELS & START VALUES as objects(to ease specification)
LabThMZf <-c('Tmz_11f','imz_11f')	# THs for var 1 and 2 for a twin individual (mz)
LabThMZm <-c('Tmz_11m','imz_11m')
LabCorMZf	<-c('r21f','rMZ1f','MZxtxtf','MZxtxtf','rMZ2f','r21f') #labels for correlation matrix
LabCorMZm <-c('r21m','rMZ1m','MZxtxtm','MZxtxtm','rMZ2m','r21m')

LabThDZf	<-c('Tdz_11f','idz_11f')
LabThDZm <-c('Tdz_11m','idz_11m')
LabCorDZf	<-c('r21f','rDZ1f','DZxtxtf','DZxtxtf','rDZ2f','r21f')
LabCorDZm <-c('r21m','rDZ1m','DZxtxtm','DZxtxtm','rDZ2m','r21m')

LabCovf	<-c('Thnnsif', 'BageThnnsif', 'BageThmjf', 'BageThmjf')
LabCovm <-c('BageThnnsim', 'BageThnnsim', 'BageTmjm', 'BageTmjm')

ThPat <-c(T,T)

# StCorMZf	<-c(.1, .6, .2, .2, .6, .1)
# StCorMZm <-c(.1, .6, .2, .2, .6, .1)
# StCorDZf	<-c(.1, .2, .1, .1, .2, .1)
# StCorDZm <-c(.1, .2, .1, .1, .2, .1)
StCorMZf  <-c(rp, h1sq, covG, covG, h2sq, rp)
StCorMZm <-c(rp, h1sq, covG, covG, h2sq, rp)
StCorDZf	<-c(rp, .5*h1sq, .5*covG, .5*covG, .5*h2sq, rp)
StCorDZm <-c(rp, .5*h1sq, .5*covG, .5*covG, .5*h2sq, rp)

#StTH	<-c(1.2816, 0.9741)
StTH <- thresh[1:2]

# Define definition variables to hold the Covariates
obsAge1  <- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age1"), name="Age1")
obsAge2	<- mxMatrix( type="Full", nrow=1, ncol=1, free=F, labels=c("data.age2"), name="Age2")

# Matrix & Algebra for expected means (SND), Thresholds, effect of Age on Th and correlations
Mean	<-mxMatrix( type="Zero", nrow=1, ncol=ntv, name="M" )

betaAf	<-mxMatrix( type="Full", nrow=nth, ncol=nv, free=T, values=.2, labels=LabCovf, name="BageTHf" )
#betaAm <-mxMatrix( type="Full", nrow=nth, ncol=nv, free=T, values=.2, labels=LabCovm, name="BageTHm" )
betaAm <-mxMatrix( type="Full", nrow=nth, ncol=nv, free=T, values=.2, labels=c("Thnnsim", "BageThnnsim"), name="BageTHm" )

inc <-mxMatrix( type="Lower",nrow=nth, ncol=nth, free=F, values=1, name="Low")

Tmzf	<-mxMatrix( type="Full", nrow=nth, ncol=nv, free=ThPat, values=StTH, lbound= -3, ubound=3, labels=LabThMZf, name="ThMZf")
Tmzm <-mxMatrix( type="Full", nrow=nth, ncol=nv, free=ThPat, values=StTH, lbound= -3, ubound=3, labels=LabThMZm, name="ThMZm") #BMH now give only one lower limit
ThresMZf	<-mxAlgebra( expression= cbind(Low%*%ThMZf + BageTHf%x%Age1, Low%*%ThMZf + BageTHf%x%Age2), name="expThresMZf")
ThresMZm <-mxAlgebra( expression= cbind(Low%*%ThMZm + BageTHm%x%Age1, Low%*%ThMZm + BageTHm%x%Age2), name="expThresMZm")
CorMZf	<-mxMatrix( type="Stand", nrow=ntv, ncol=ntv, free=T, values=StCorMZf, labels=LabCorMZf, lbound=-.99, ubound=.99, name="expCorMZf")
CorMZm <-mxMatrix( type="Stand", nrow=ntv, ncol=ntv, free=T, values=StCorMZm, labels=LabCorMZm, lbound=-.99, ubound=.99, name="expCorMZm")

Tdzf	<-mxMatrix( type="Full", nrow=nth, ncol=nv, free=ThPat, values=StTH, lbound= -3, ubound=3, labels=LabThDZf, name="ThDZf") #BMH now give only one lower limit
Tdzm <-mxMatrix( type="Full", nrow=nth, ncol=nv, free=ThPat, values=StTH, lbound= -3, ubound=3, labels=LabThDZm, name="ThDZm")
ThresDZf	<-mxAlgebra( expression= cbind(Low%*%ThDZf + BageTHf%x%Age1, Low%*%ThDZf + BageTHf%x%Age2), name="expThresDZf")
ThresDZm <-mxAlgebra( expression= cbind(Low%*%ThDZm + BageTHm%x%Age1, Low%*%ThDZm + BageTHm%x%Age2), name="expThresDZm")
CorDZf	<-mxMatrix( type="Stand", nrow=ntv, ncol=ntv, free=T, values=StCorDZf, labels=LabCorDZf, lbound=-.99, ubound=.99, name="expCorDZf")
CorDZm <-mxMatrix( type="Stand", nrow=ntv, ncol=ntv, free=T, values=StCorDZm, labels=LabCorDZm, lbound=-.99, ubound=.99, name="expCorDZm")

# Data objects for Multiple Groups
dataMZf <- mxData( observed=mzfData, type="raw" )
dataMZm <- mxData( observed=mzmData, type="raw" )
dataDZf <- mxData( observed=dzfData, type="raw" )
dataDZm <- mxData( observed=dzmData, type="raw" )

expMZf <- mxExpectationNormal( covariance="expCorMZf", means="M", dimnames=selVars, thresholds="expThresMZf" )
expMZm <- mxExpectationNormal( covariance="expCorMZm", means="M", dimnames=selVars, thresholds="expThresMZm" )
expDZf <- mxExpectationNormal( covariance="expCorDZf", means="M", dimnames=selVars, thresholds="expThresDZf" )
expDZm <- mxExpectationNormal( covariance="expCorDZm", means="M", dimnames=selVars, thresholds="expThresDZm" )

#FitFunctions
fun <- mxFitFunctionML()
fitFunction <- mxFitFunctionMultigroup(c("MZf.fitfunction", "MZm.fitfunction", "DZf.fitfunction", "DZm.fitfunction"))

# Combine Groups
modelMZf	<- mxModel( obsAge1, obsAge2, Mean, betaAf, Tmzf, inc, ThresMZf, CorMZf, dataMZf, expMZf, fun, name="MZf" )
modelMZm <- mxModel( obsAge1, obsAge2, Mean, betaAm, Tmzm, inc, ThresMZm, CorMZm, dataMZm, expMZm, fun, name="MZm" )
modelDZf	<- mxModel( obsAge1, obsAge2, Mean, betaAf, Tdzf, inc, ThresDZf, CorDZf, dataDZf, expDZf, fun, name="DZf" )
modelDZm <- mxModel( obsAge1, obsAge2, Mean, betaAm, Tdzm, inc, ThresDZm, CorDZm, dataDZm, expDZm, fun, name="DZm" )

SatModel <- mxModel( "Sat", modelMZf, modelMZm, modelDZf, modelDZm, fitFunction)

# my code: where the madness starts
SatModel$MZf$Age1
SatModel$MZf$BageTHf
mxEval(MZf.BageTHf %x% MZf.Age1, SatModel, compute=T)
mxEval(matrix(1, 1, 1) %*% MZf.Age1, SatModel, compute=T)

# -------------------------------------------------------------------------------------------------------------------------------
# 1) RUN Saturated Model

SatFit	<- mxRun(SatModel, intervals=F)
summary(SatFit)

# more of my code
mod1 <- mxModel(obsAge1, obsAge2, Mean, betaAf, Tmzf, inc,
                ThresMZf, CorMZf, dataMZf, expMZf, fun,
                name="MZf", mxFitFunctionML())
mod1Res <- mxRun(mod1)

# definition variables and dependencies should be set to NA upon return
omxCheckEquals(mod1Res$Age1$values, dataMZf$observed[1,'age1'])
omxCheckEquals(mod1Res$Age2$values, dataMZf$observed[1,'age2'])
omxCheckCloseEnough(c(mod1Res$ThMZf$values), c(0.732, 0.305), .01)
omxCheckCloseEnough(c(mod1Res$BageTHf$values), c(-0.008, 0), .01)
omxCheckEquals(nrow(mod1Res$data$observed), 250)
omxCheckEquals(mod1Res$expThresMZf$result,
               mxEval(expThresMZf, mod1Res, compute=T, defvar.row = 1))
