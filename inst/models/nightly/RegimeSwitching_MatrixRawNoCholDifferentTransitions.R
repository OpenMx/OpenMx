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

# -----------------------------------------------------------------------------
# Program: RegimeSwitching_MatrixRaw.R  
# Author: Michael Neale
# Date: 2010.09.17 
#
# ModelType: Growth Mixture with Switching
# DataType: Continuous
# Field: None
#
# Purpose: 
#       Regime Switching Growth Mixture Model
#       Matrix style model input - Raw data input
#       Borrows from seminal work by Conor Dolan et al 2005 in SEM Journal
#
# RevisionHistory:
# -----------------------------------------------------------------------------

# Load Libraries
require(OpenMx)
#mxOption(NULL, "Number of Threads", 8L)
#mxOption(NULL, "Default optimizer", "NPSOL")
set.seed(1)

"%&%" <- OpenMx::"%&%"  # ensure we don't use the %&% from Matrix

# -----------------------------------------------------------------------------

readData <- function(path) {
	read.table(file=path,header=FALSE,na.strings='-9.00',col.names=c(paste("time",1:4,sep="")))
}

# Prepare Data
dolan2005Data <- suppressWarnings(try(readData('data/sel.txt')))
if (is(dolan2005Data, "try-error")) dolan2005Data <- readData('models/nightly/data/sel.txt')
# -----------------------------------------------------------------------------

# Number of occasions, factors & regimes
nocc <- 4
nfac <- 2
nregime <- 3
nfacxnregime <- nfac * nregime
# -----------------------------------------------------------------------------

# Set up labels for the models; these populate larger vectors
baseMeanLabels<-c("meanIntercept","meanSlope")
baseCovarianceLabels<-c("varIntercept","covInterceptSlope","covInterceptSlope","varSlope")
# -----------------------------------------------------------------------------

# Set up symmetric factor covariance matrix - same thing for all nregime^nocc models
factorCovariances <- mxMatrix(type="Symm",nfacxnregime,nfacxnregime,name="factorCovariances")

j <- 1
for (i in 1:nregime)
{ 
    factorCovariances$labels[j:(j+1),j:(j+1)] <- paste(baseCovarianceLabels,i,sep="")
    j <- j + nfac
}
factorCovariances$free[2,2] <- TRUE
factorCovariances$free[6,6] <- TRUE
factorCovariances$lbound[2,2] <- .0001
factorCovariances$lbound[6,6] <- .0001
factorCovariances$values[2,2] <- .1
factorCovariances$values[6,6] <- .05
#
# -----------------------------------------------------------------------------

# Set up base factor loading matrix; its values change with each regime-switched models
factorLoadings <- mxMatrix(type="Full", nrow=nocc, ncol=nfacxnregime, values=c(rep(1,nocc),0:(nocc-1),rep(0,(nocc*(nfacxnregime-2)))), free=F, name="factorLoadings")
# -----------------------------------------------------------------------------

# Factor Mean matrix; same for all models.  Funky method of getting cartesian product of baseLabels and 1:nregime
factorMeans <- mxMatrix(type="Full", nrow=nfacxnregime, ncol=1, free=c(T,T,T,F,T,T), values=c(3.35,0.45,10.1,0,1.17,0.22), name="factorMeans")
factorMeans$labels[] <- apply(expand.grid(baseMeanLabels,1:nregime),1,function(x) paste(x,collapse=""))
# -----------------------------------------------------------------------------

# Set up base residual variance matrix; its parameter labels change with each regime-switched model
residuals <- mxMatrix(type="Diag", nrow=nocc, ncol=nocc, free=T, values=1, lbound=.1, labels="residual1", name="residuals")
# -----------------------------------------------------------------------------


# Create a generic MxModel object for Regime 1, allowing for nregime * nfac latent variables which influence nocc observed variables
# It is done this way to allow for possible covariances between Regime latent growth curve factors
regime1 <- mxModel("regime1", factorCovariances, residuals, factorMeans, factorLoadings,
    mxAlgebra(expression=factorLoadings %&% (factorCovariances) + residuals, name="expectedCovariances"),
    mxAlgebra(expression=t(factorLoadings %*% factorMeans), name="expectedMeans"),
    mxExpectationNormal(covariance="expectedCovariances", means="expectedMeans", dimnames = names(dolan2005Data)),
		   mxFitFunctionML(vector=TRUE)
                 )
# -----------------------------------------------------------------------------

# Create a list of MxModel objects and their names for regimeswitches 1...nregime^nocc by overwriting the parameter labels 
#   of the factor means and covariances of class 1's model
modelList <- vector("list",(nregime^nocc))
modelNames <- vector("list",(nregime^nocc))
# -----------------------------------------------------------------------------

# Construct a nregime^nocc by nocc matrix with the relevant switching patterns in it
switchPatterns<-matrix(0,nregime^nocc,0)
for (j in 1:nocc) {
    conor <- matrix(1,nregime,nocc)
    conor[,j]<-1:nregime
    tmp1 <- conor[,1]
    for (k in 2:nocc) {
        tmp1 <- tmp1 %x% conor[,k]
    }
    switchPatterns<-cbind(switchPatterns,tmp1)
    }
# -----------------------------------------------------------------------------

# Use the switching patterns to modify the models in the list
# First set up a couple matrices to hold parameter label numbers and starting values
resNumber <- matrix(c(rep((1:nregime),nocc)),nrow=nregime,ncol=nocc)
resValues <- matrix(rep(c(2,5,3)),nrow=nregime,ncol=nocc)
# Then another one for the basis functions (factor loadings)
ncomponents <- 2
# NB, lots of things need to change, including line below, if ncomponents >2, particularly factorCovariances and factorLoadings
basisFunctions <- matrix( c(rep(1,nocc),0:(nocc-1)),nocc,ncomponents)
# This loop makes use of a patternVector z to describe the particular combination of regimes that a person may be in
ii <- 0
for (i in 1:(nregime^nocc))
{
    temp <- regime1
    patternVector <- t(switchPatterns[i,])
    for (j in 1:nocc) 
        {
        ii <- ii+1
        z <- matrix(0,1,nregime)
        z[patternVector[j]] <- 1
        temp$residuals$labels[j,j] <- paste("residual",resNumber[patternVector[j],j],sep="")
        temp$residuals$values[j,j] <- resValues[patternVector[j],j]
        temp$factorLoadings$values[j,] <- z %x% basisFunctions[j,]
        }

    temp <- mxModel(name=paste("Regime",i,sep=""), temp)
    modelNames[i] <- paste("Regime",i,sep="")
    modelList[i] <- list(temp)
}
# -----------------------------------------------------------------------------
# Now we construct the weight matrix algebra
# in Classic Mx it was:
# g unit nr 1              	
# p full nr 1 fr				! probs initial
# q full nr nr fr    			! transition probs 
# end matrices;
# !
# begin algebra;              	
# d=\m2v(q');    			!  transition probs into vector
# ! mixing proportions 
# ! a=(p$g$g).(d$g).(g$d)					! nt=3
# a=(p$g$g$g).(d$g$g).(g$d$g).(g$g$d);		! nt=4
# and I am not going to f with making it general just now...

# To avoid constraints, we express initial and transition matrices as proportions by bounding them to be positive and expressing as proportions (TBD)
rawInitialProbs <- mxMatrix(type="Full", nrow=nregime, ncol=1, free=c(F,rep(T,(nregime-1))), values=c(.35,.15,.50), lbound=.00001, name="rawInitialProbs")
rawTransitionProbs1 <- mxMatrix(type="Full", nrow=nregime, ncol=nregime, free=c(F,rep(T,(nregime-1))), values=c(1,.2,.1,1,.7,.05,1,.05,.65), lbound=.00001, name="rawTransitionProbs1")
initialProbs <- mxAlgebra(rawInitialProbs %x% (1/sum(rawInitialProbs)), name="initialProbs")
# was
#transitionProbs1 <- mxAlgebra(rawTransitionProbs1 %x% (1/sum(rawTransitionProbs1)), name="transitionProbs1")
transitionProbs1 <- mxAlgebra(rawTransitionProbs1 / (unitNregime %x% (t(unitNregime) %*% rawTransitionProbs1)), name="transitionProbs1")
transitionVector1 <- mxAlgebra(cvectorize(transitionProbs1), name="transitionVector1")

rawTransitionProbs2 <- mxMatrix(type="Full", nrow=nregime, ncol=nregime, free=c(F,rep(T,(nregime-1))), values=c(1,.2,.1,1,.7,.05,1,.05,.65), lbound=.00001, name="rawTransitionProbs2")
# was
#transitionProbs2 <- mxAlgebra(rawTransitionProbs2 %x% (1/sum(rawTransitionProbs2)), name="transitionProbs2")
transitionProbs2 <- mxAlgebra(rawTransitionProbs2 / (unitNregime %x% (t(unitNregime) %*% rawTransitionProbs2)), name="transitionProbs2")
transitionVector2 <- mxAlgebra(cvectorize(transitionProbs2), name="transitionVector2")

rawTransitionProbs3 <- mxMatrix(type="Full", nrow=nregime, ncol=nregime, free=c(F,rep(T,(nregime-1))), values=c(1,.2,.1,1,.7,.05,1,.05,.65), lbound=.00001, name="rawTransitionProbs3")
# was
#transitionProbs3 <- mxAlgebra(rawTransitionProbs3 %x% (1/sum(rawTransitionProbs3)), name="transitionProbs3")
transitionProbs3 <- mxAlgebra(rawTransitionProbs3 / (unitNregime %x% (t(unitNregime) %*% rawTransitionProbs3)), name="transitionProbs3")
transitionVector3 <- mxAlgebra(cvectorize(transitionProbs3), name="transitionVector3")


unitNregime <- mxMatrix(type="Unit",nrow=nregime, ncol=1, name="unitNregime")
weights <- mxAlgebra((initialProbs %x% unitNregime %x% unitNregime %x% unitNregime) * (transitionVector1 %x% unitNregime %x% unitNregime) * (unitNregime %x% transitionVector2 %x% unitNregime) * (unitNregime %x% unitNregime %x% transitionVector3), name="weights")
# -----------------------------------------------------------------------------

# Build up objective function using generous doses of paste
objectives <- paste(modelNames, 'objective', sep=".")
modelnumbers <- 1:(nregime^nocc)

components <- paste("weights[", modelnumbers, ",1]", " %x% ", objectives, sep = '')
componentsSum <- paste(components, collapse = " + ")
algebraString <- paste("mxAlgebra(-2*sum(log(", componentsSum, ")), name='mixtureObj')", sep = "")
algebraObjective <- eval(parse(text=algebraString)[[1]])
obj <- mxFitFunctionAlgebra("mixtureObj")
# -----------------------------------------------------------------------------
      
# Construct MxModel from all the various pieces; add options and then run it and summarize results
rsgcmmDiffT <- mxModel("Regime Switching Growth Curve Model", mxData(observed=dolan2005Data, type="raw"), modelList, rawInitialProbs, initialProbs, rawTransitionProbs1, transitionProbs1, transitionVector1, rawTransitionProbs2, transitionProbs2, transitionVector2, rawTransitionProbs3, transitionProbs3, transitionVector3, unitNregime, weights, algebraObjective, obj)      
rsgcmmDiffT <- mxOption(rsgcmmDiffT, 'Checkpoint Count', 1)
#rsgcmm <- mxOption(rsgcmm, "Standard Errors", "No")
#rsgcmm <- mxOption(rsgcmm, "Calculate Hessian", "No")

rsgcmmDiffTFit <- mxTryHard(rsgcmmDiffT, extraTries=4)

summary(rsgcmmDiffTFit)

# CSOLNP 10095
# SLSQP 9933.171 
# NPSOL 9962.865 

# -----------------------------------------------------------------------------

# Check results to see if they are within specified bounds
omxCheckCloseEnough(rsgcmmDiffTFit$output$Minus2LogLikelihood, 9933, 170)
#omxCheckCloseEnough(max(mxEval(rsgcmmFit$output$transitionProbs, rsgcmmFit)), 0.4790, 0.01)
#omxCheckCloseEnough(min(mxEval(rsgcmmFit$output$transitionProbs, rsgcmmFit)), 0.1608, 0.01)
## -----------------------------------------------------------------------------
