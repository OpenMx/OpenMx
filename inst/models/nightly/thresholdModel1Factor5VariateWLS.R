#
#   Copyright 2007-2018 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


#
#  OpenMx Ordinal Data Example
#  Revision history:
#		Michael Neale 14 Aug 2010
#

# Step 1: load libraries
require(OpenMx)

#
# Step 2: set up simulation parameters 
# Note: nVariables>=3, nThresholds>=1, nSubjects>=nVariables*nThresholds (maybe more)
# and model should be identified
#
nVariables<-5
nFactors<-1
nThresholds<-3
nSubjects<-500
isIdentified<-function(nVariables,nFactors) as.logical(1+sign((nVariables*(nVariables-1)/2) -  nVariables*nFactors + nFactors*(nFactors-1)/2))
# if this function returns FALSE then model is not identified, otherwise it is.
isIdentified(nVariables,nFactors)

loadings <- matrix(.7,nrow=nVariables,ncol=nFactors)
residuals <- 1 - (loadings * loadings)
sigma <- loadings %*% t(loadings) + vec2diag(residuals)
mu <- matrix(0,nrow=nVariables,ncol=1)
# Step 3: simulate multivariate normal data
set.seed(1234)
continuousData <- mvtnorm::rmvnorm(n=nSubjects,mu,sigma)

# Step 4: chop continuous variables into ordinal data 
# with nThresholds+1 approximately equal categories, based on 1st variable
quants<-quantile(continuousData[,1],  probs = c((1:nThresholds)/(nThresholds+1)))
ordinalData<-matrix(0,nrow=nSubjects,ncol=nVariables)
for(i in 1:nVariables)
{
ordinalData[,i] <- cut(as.vector(continuousData[,i]),c(-Inf,quants,Inf))
}

# Step 5: make the ordinal variables into R factors
ordinalData <- mxFactor(as.data.frame(ordinalData),levels=c(1:(nThresholds+1)))

# Step 6: name the variables
fruitynames<-paste("banana",1:nVariables,sep="")
names(ordinalData)<-fruitynames


thresholdModel <- mxModel("thresholdModel",
	mxMatrix("Full", nVariables, nFactors, values=0.2, free=TRUE, lbound=-.99, ubound=.99, name="L"),
	mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
	mxMatrix("Zero", 1, nVariables, name="M"),
	mxAlgebra(vectorofOnes - (diag2vec(L %*% t(L))) , name="E"),
	mxAlgebra(L %*% t(L) + vec2diag(E), name="impliedCovs"),
	mxMatrix("Full", 
            name="thresholdDeviations", nrow=nThresholds, ncol=nVariables,
            values=.2,
            free = TRUE, 
            lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , nVariables),
            dimnames = list(c(), fruitynames)),
    mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
    mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
            mxFitFunctionML(),mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = fruitynames, thresholds="thresholdMatrix"),
            mxData(observed=ordinalData, type='raw')
)

thresholdModelrun <- mxRun(thresholdModel)
thresholdSaturated <- mxRefModels(thresholdModelrun, run=TRUE)
summary(thresholdModelrun, refModels=thresholdSaturated)

a <- Sys.time()
thresholdModelAuto <- mxAutoStart(thresholdModel)
b <- Sys.time()
b-a #about 2 seconds on my laptop

thresholdModelAutoRun <- mxRun(thresholdModelAuto)

(b-a) + summary(thresholdModelAutoRun)$wallTime
summary(thresholdModelrun)$wallTime
# 37 sec for the automatic start values time plus the estimation time from the auto starts
# 64 sec for the estimation time from the user starts
# auto starts provides a net boost in performance

a <- proc.time()
thresholdModelWLS <- mxModel(thresholdModel, name="WLSThresholdModel", mxDataWLS(ordinalData, type="ULS"), #Change type here!!!
	mxMatrix('Zero', nrow=1, ncol=nVariables, name='impliedMeans'),
	mxExpectationNormal(covariance="impliedCovs", means='impliedMeans', dimnames = fruitynames, thresholds="thresholdMatrix"),
	mxFitFunctionWLS())
thresholdModelWLSrun <- mxRun(thresholdModelWLS)
b <- proc.time()
b-a
summary(thresholdModelrun)$wallTime

summary(thresholdModelWLSrun)

wls.L <- mxEval(L, thresholdModelWLSrun) #should be all 0.7
wls.T <- mxEval(thresholdMatrix, thresholdModelWLSrun) #should be all quants

ml.L <- mxEval(L, thresholdModelrun) #should be all 0.7
ml.T <- mxEval(thresholdMatrix, thresholdModelrun) #should be all quants

auto.L <- mxEval(L, thresholdModelAuto)
auto.T <- mxEval(thresholdMatrix, thresholdModelAuto, compute=TRUE)

rms <- function(x, y){sqrt(mean((x-y)^2))}

omxCheckTrue(rms(wls.L, .7) < 0.025)
rms(ml.L, .7)
omxCheckTrue(rms(ml.L, auto.L) < 0.01)

omxCheckTrue(rms(wls.T, quants) < 0.05)
rms(ml.T, quants)

omxCheckTrue(rms(wls.L, auto.L) < 1e-6)
omxCheckTrue(rms(wls.T, auto.T) < 1e-6)

ml.sum <- summary(thresholdModelrun, refModels=thresholdSaturated)
wls.sum <- summary(thresholdModelWLSrun)
omxCheckWithinPercentError(wls.sum$Chi, 0.653, percent=10)
omxCheckWithinPercentError(ml.sum$Chi, wls.sum$Chi, percent=16)
omxCheckEquals(ml.sum$ChiDoF, wls.sum$ChiDoF)

ciModel <- mxModel(thresholdModelWLSrun, mxCI("L"))
omxCheckError(mxRun(ciModel, intervals=TRUE), "Confidence intervals are not supported for units r'Wr")


#------------------------------------------------------------------------------
tmod2 <- mxModel("thresholdModel2",
	mxMatrix("Full", nVariables, nFactors, values=0.2, free=TRUE, lbound=-.99, ubound=1.5, name="L"),
	mxMatrix("Diag", nVariables, nVariables, values=.1, free=TRUE, lbound=1e-8, name="R"),
	mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
	mxMatrix("Full", 1, nVariables, values=0, free=TRUE, name="M"),
	mxAlgebra(L %*% t(L) + R, name="impliedCovs"),
	mxMatrix("Full", nThresholds, nVariables, values=c(0, 1, 2), name="Thresh"),
	mxFitFunctionML(),
	mxExpectationNormal(covariance="impliedCovs", means="M", dimnames = fruitynames, thresholds="Thresh"),
	mxData(observed=ordinalData, type='raw')
)

trun2 <- mxRun(tmod2)

a <- proc.time()
wmod2 <- mxModel(tmod2, mxDataWLS(ordinalData), mxFitFunctionWLS(),
	mxAlgebra(cov2cor(impliedCovs), name='newCov'),
	mxMatrix("Unit", nrow=nThresholds, ncol=1, name="UnitVector"),
	mxAlgebra(UnitVector %x% t(sqrt(diag2vec(impliedCovs))), name='theStandardDeviations'),
	mxAlgebra(UnitVector %x% M, name='theM'),
	mxAlgebra( (Thresh-theM)/theStandardDeviations, name='newThresh'),
	mxMatrix('Zero', 1, nVariables, name='MZ'),
	mxExpectationNormal(covariance='newCov', means='MZ', thresholds='newThresh', dimnames = fruitynames) #N.B. means left out on purpose
)


#mxEval(theM, wmod2, compute=TRUE)
#mxEval(Thresh, wmod2, compute=TRUE)
#mxEval(theStandardDeviations, wmod2, compute=TRUE)
#mxEval(newCov, wmod2, compute=TRUE)
#mxEval(newThresh, wmod2, compute=TRUE)

wrun2 <- mxRun(wmod2)
b <- proc.time()

b-a
summary(trun2)$wallTime

cbind(omxGetParameters(trun2), omxGetParameters(wrun2))

plot(omxGetParameters(trun2), omxGetParameters(wrun2))
abline(a=0, b=1)

omxCheckCloseEnough(rms(omxGetParameters(trun2), omxGetParameters(wrun2)), 0, .035)
omxCheckCloseEnough(cor(omxGetParameters(trun2), omxGetParameters(wrun2)), 1, .06)


# new style for model 2
wmod2a <- mxModel(tmod2, mxDataWLS(ordinalData), mxFitFunctionWLS())
wrun2a <- mxRun(wmod2a)

cbind(omxGetParameters(trun2), omxGetParameters(wrun2), omxGetParameters(wrun2a))

# Check that old/hard stadardization and new/easy standardization give the same
#  answer.
omxCheckCloseEnough(omxGetParameters(wrun2), omxGetParameters(wrun2a), 1e-4)

#------------------------------------------------------------------------------

require(OpenMx)
a <- factor(sample(c('a', 'b', 'c'), size=100, replace=T))
b <- factor(a, levels=c('a', 'b', 'c', 'd')) #create factor is unused level 'd'
ma <- mxFactor(a, levels=levels(a))
mb <- mxFactor(b, levels=levels(b))
#mb[mb %in% 'c'] <- 'd' #make 'c' the unused level instead of 'd'
ds <- data.frame(a=ma, b=mb)

omxCheckError(wd <- mxDataWLS(ds), "Variable 'b' has a zero frequency category 'd'.
Eliminate this level in your mxFactor() or combine categories in some other way.
Do not pass go. Do not collect $200.")

