#
#   Copyright 2007-2016 The OpenMx Project
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


#------------------------------------------------------------------------------
# Frontpage model in LISREL form
#  with factor scores

require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("OneFactor", 
                       type="LISREL",
                       manifestVars=list(exo=manifests), 
                       latentVars=list(exo=latents),
                       mxPath(from=latents, to=manifests),
                       mxPath(from=manifests, arrows=2),
                       mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
                       mxPath(from='one', to=manifests),
                       mxData(observed=demoOneFactor, type="raw"))
summary(factorRun <- mxRun(factorModel))


r1 <- mxFactorScores(factorRun, 'ML')
r2 <- mxFactorScores(factorRun, 'Regression')
r3 <- mxFactorScores(factorRun, 'WeightedML')


omxCheckCloseEnough(cor(r1[,,1], r2[,,1]), 1)


#pdf('plotFactorScores.pdf')
#plot(r2[,,1], r1[,,1], main='Factor Scoring Methods', xlab='Regression Score', ylab='Likelihood Score')
#legend('bottomright', legend=c('ML', 'Weighted ML'), pch=1, col=c('black', 'blue'), lty=1)
#points(r2[,,1], r3[,,1], col='blue')
#abline(v=0, h=0)
lmS <- lm(r3[,,1] ~ r2[,,1])
lmR <- lm(r1[,,1] ~ r2[,,1])
#abline(lmS, col='blue')
#abline(lmR, col='black')
#dev.off()

summary(lmS)
summary(lmR)

omxCheckCloseEnough(coef(lmS), c(0, .5), 0.01)
omxCheckCloseEnough(coef(lmR), c(0, 1), 0.03)


#Test warning about no SEs:
#factorRun <- mxOption(factorRun,"Standard Errors","No")
mxOption(NULL,"Standard Errors","No")
omxCheckWarning(
	mxFactorScores(factorRun,"ML"),
	"factor-score standard errors not available from MxModel 'OneFactor' because calculating SEs is turned off for that model (possibly due to one or more MxConstraints)")
omxCheckWarning(
	mxFactorScores(factorRun,"WeightedML"),
	"factor-score standard errors not available from MxModel 'OneFactor' because calculating SEs is turned off for that model (possibly due to one or more MxConstraints)")
mxOption(NULL,"Standard Errors","Yes")



#TODO compare standard errors
#cbind(r2[,,2], r1[,,2], r3[,,2])




#------------------------------------------------------------------------------
# Ordinal example similar to thresholdModel1Factor3Variate.R

# Step 1: load libraries
require(OpenMx)

#
# Step 2: set up simulation parameters 
# Note: nVariables>=3, nThresholds>=1, nSubjects>=nVariables*nThresholds (maybe more)
# and model should be identified
#
nVariables<-3
nFactors<-1
nThresholds<-3
nSubjects<-500
isIdentified<-function(nVariables,nFactors) as.logical(1+sign((nVariables*(nVariables-1)/2) -  nVariables*nFactors + nFactors*(nFactors-1)/2))
# if this function returns FALSE then model is not identified, otherwise it is.
isIdentified(nVariables,nFactors)



# Step 3: simulate multivariate normal data
set.seed(1234)

quants <- quantile(rnorm(1001),  probs = c((1:nThresholds)/(nThresholds+1)))
invL <- matrix(0, nThresholds, nThresholds)
invL[lower.tri(invL, diag=TRUE)] <- 1
thresholdStart <- solve(invL) %*% quants


fruitynames<-paste("banana",1:nVariables,sep="")


lis <- mxModel("thresholdModel",
	mxMatrix("Full", nVariables, nFactors, values=0.7, free=TRUE, lbound=-.99, ubound=.99, name="L", dimnames=list(fruitynames, 'F')),
	mxMatrix("Unit", nVariables, 1, name="vectorofOnes"),
	mxMatrix("Zero", nVariables, 1, name="M", dimnames=list(fruitynames, 'F')),
	mxAlgebra(vec2diag(vectorofOnes - (diag2vec(L %*% t(L)))) , name="E"),
	mxMatrix('Iden', 1, 1, name='P', dimnames=list('F', 'F')),
	mxMatrix('Zero', 1, 1, name='K', dimnames=list('F', c())),
	mxMatrix("Full", 
            name="thresholdDeviations", nrow=nThresholds, ncol=nVariables,
            values=thresholdStart,
            free = TRUE, 
            lbound = rep( c(-Inf,rep(.01,(nThresholds-1))) , nVariables),
            dimnames = list(c(), fruitynames)),
    mxMatrix("Lower",nThresholds,nThresholds,values=1,free=F,name="unitLower"),
    mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
            mxFitFunctionML(),mxExpectationLISREL(LX='L', TX='M', PH='P', TD='E', KA='K', thresholds="thresholdMatrix")
)

# Generate data based on model
ordinalData <- mxGenerateData(lis, nSubjects)

# Put the simulated data in the model and run
lis <- mxModel(lis, mxData(observed=ordinalData, type='raw'))
lisr <- mxRun(lis)


# Compute factor scores for the model
lism <- mxFactorScores(lisr)
lisw <- mxFactorScores(lisr, 'WeightedML')

omxCheckError(lisreg <- mxFactorScores(lisr, 'Regression'), "Regression factor scores cannot be computed when there are thresholds (ordinal data).")

omxCheckCloseEnough(cor(lism[,,1], lisw[,,1]), .868, 0.01)

#pdf('plotOrdinalFactorScores.pdf')
#plot(lism[,,1], main='Ordinal Factor Scoring Methods', xlab='Sorted Data Row', ylab='Factor Score')
#points(lisw[,,1], col='blue')
#legend('bottomright', legend=c('ML', 'Weighted ML', paste('r =', round(cor(lism[,,1], lisw[,,1]), 2))), pch=c(1, 1, NA), col=c('black', 'blue', 'black'))
#dev.off()


#------------------------------------------------------------------------------
# Frontpage model in RAM form

require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("OneFactor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxPath(from='one', to=manifests),
      mxData(observed=demoOneFactor, type="raw"))
summary(factorRamRun <- mxRun(factorModel))

rr1 <- mxFactorScores(factorRamRun, 'ML')
rr2 <- mxFactorScores(factorRamRun, 'Regression')
rr3 <- mxFactorScores(factorRamRun, 'WeightedML')


omxCheckCloseEnough(cor(rr1[,,1], r1[,,1]), 1)
omxCheckCloseEnough(cor(rr2[,,1], r2[,,1]), 1)
omxCheckCloseEnough(cor(rr3[,,1], r3[,,1]), 1)

rms <- function(x, y){sqrt(mean((x-y)^2))}

omxCheckCloseEnough(rms(rr1[,,1], r1[,,1]), 0, .1)
omxCheckCloseEnough(rms(rr2[,,1], r2[,,1]), 0, .1)
omxCheckCloseEnough(rms(rr3[,,1], r3[,,1]), 0, .1)


#------------------------------------------------------------------------------
# Multiple group example

factorModel2 <- mxModel(factorModel, name='ConstrainedResidual',
	mxPath(from=manifests, arrows=2, labels='resid')
)

twoGroup <- mxModel(model='multipleGroup', factorModel, factorModel2,
	mxFitFunctionMultigroup(c('OneFactor', 'ConstrainedResidual'))
)
twoGroupRun <- mxRun(twoGroup)

tg1 <- mxFactorScores(twoGroupRun, 'ML')
tg2 <- mxFactorScores(model=twoGroupRun,"Regression")

omxCheckCloseEnough(rms(tg1[[1]], rr1), 0, .001)
omxCheckCloseEnough(rms(tg1[[2]], rr1), 0, .01)
omxCheckCloseEnough(rms(tg2[[1]], rr2), 0, .001)
omxCheckCloseEnough(rms(tg2[[2]], rr2), 0, .01)


#-------------------------------------------------------------------
#Ensure regression RAM scoring does not fail in the presence of missing data:
demoOneFactor.miss <- as.matrix(demoOneFactor)
demoOneFactor.miss[sample(1:2500,size=100,replace=F)] <- NA
demoOneFactor.miss[100,] <- NA
factorModel.miss <- mxModel("OneFactor",
											 type="RAM",
											 manifestVars = manifests,
											 latentVars = latents,
											 mxPath(from=latents, to=manifests),
											 mxPath(from=manifests, arrows=2),
											 mxPath(from=latents, arrows=2,
											 			 free=FALSE, values=1.0),
											 mxPath(from='one', to=manifests),
											 mxData(observed=demoOneFactor.miss, type="raw"))
factorRamRun.miss <- mxRun(factorModel.miss)
regs <- mxFactorScores(model=factorRamRun.miss,"Regression")
omxCheckEquals(regs[100,1,1],0)
omxCheckEquals(regs[100,1,2],1)
omxCheckTrue(cor(regs[,,1], rr2[,,1]) > 0.95)

