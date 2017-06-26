#------------------------------------------------------------------------------
# Load required packages

require(OpenMx)

#------------------------------------------------------------------------------
# Read-in and list-ify data

data(myLongitudinalData)
matplot(t(myLongitudinalData), type='l')
dataL <- list()
for(i in 1:nrow(myLongitudinalData)){
	obsI <- unlist(myLongitudinalData[i,])
	dataL[[i]] <- data.frame(y=obsI, ID=i, time=0:4)
}
nSubj <- length(dataL)


#------------------------------------------------------------------------------
# Create state space matrices for latent growth

amat <- mxMatrix("Iden", 2, 2, name="A")
bmat <- mxMatrix("Zero", 2, 1, name="B")
clab <- c(NA, "data.time")
cdim <- list(c("y"), c("I", "S"))
cmat <- mxMatrix("Full", 1, 2, FALSE, c(1, NA), name="C",
    dimnames=cdim, labels=clab)
dmat <- mxMatrix("Zero", 1, 1, name="D")
qmat <- mxMatrix("Zero", 2, 2, name="Q")
rlab <- "resid"
rmat <- mxMatrix("Diag", 1, 1, TRUE, .2, name="R", labels=rlab)
xlab <- c("meanI", "meanS")
xmat <- mxMatrix("Full", 2, 1, TRUE, c(1, 1), name="x0", labels=xlab)
pval <- c(1, .5, 1)
plab <- c("varI", "covIS", "varS")
pmat <- mxMatrix("Symm", 2, 2, TRUE, pval, name="P0", plab, lbound=c(0, NA, 0))
umat <- mxMatrix("Zero", 1, 1, name="u")

modL <- list(amat, bmat, cmat, dmat, qmat, rmat, xmat,
	pmat, umat)
modNames <- paste0("Subject", 1:nSubj, "LatentGrowth")
expSS <- mxExpectationStateSpace(A="A", B="B", C="C", D="D",
	Q="Q", R="R", x0="x0", P0="P0", u="u")

#------------------------------------------------------------------------------
# Create/estimate multisubject model

indivmodels <- list()
for(k in 1:nSubj){
	DataSetForSubjectK <- dataL[[k]]
	indivmodels[[k]] <- mxModel(name=modNames[k],
		modL, expSS,
		mxFitFunctionML(),
		mxData(DataSetForSubjectK, type='raw')) 
}
multiSubjGrowth <- mxModel(name="MultiGrowth", indivmodels,
	mxFitFunctionMultigroup(modNames))

multiSubjGrowthRun <- mxRun(multiSubjGrowth)


#------------------------------------------------------------------------------
# Examine results

summary(multiSubjGrowthRun)

# Kalman Scores for Subject 1
ks <- mxKalmanScores(multiSubjGrowthRun$Subject1LatentGrowth, frontend=FALSE)

# These are equal to the empirical Bayes random effects estimates
#  from the nlme and lme4 packages
ks$xSmoothed


ksAll <- matrix(0, nrow=length(multiSubjGrowthRun$submodels), ncol=2)
for(i in 1:length(multiSubjGrowthRun$submodels)){
	ksAll[i,] <- mxKalmanScores(multiSubjGrowthRun$submodels[[i]], frontend=TRUE)$xSmoothed[1,]
}

# TODO Bug fix: when frontend=FALSE the results are different
#  and incorrect


#------------------------------------------------------------------------------
# Run the linear mixed effects model for comparison

require(lme4)
tallData <- do.call(rbind, dataL)
g <- lmer(y ~ 1 + time + (1 + time | ID), data=tallData)

gfix <- fixef(g)
gcov <- as.data.frame(VarCorr(g))$vcov
cmp <- cbind(StateSpace=coef(multiSubjGrowthRun), mixed=c(gcov[4], gfix, gcov[c(1,3,2)]))

omxCheckCloseEnough(cmp[,1], cmp[,2], 0.015)

omxCheckCloseEnough(ksAll, t(t(ranef(g)$ID) + fixef(g)), 0.01)


#------------------------------------------------------------------------------
# Create definition variable test variation with a person-level covariate
#  predicting the random means

set.seed(99)
z <- ksAll[,1] + rnorm(nrow(ksAll), mean=3.3, sd=1.8)
# implies that i = 3.3 + .5*u
lmc <- lm(ksAll[,1] ~ z)

zmat <- mxMatrix(name='Z', type='Full', nrow=2, ncol=1, labels=c(NA, 'data.z'), values=c(1, NA), free=FALSE)
emat <- mxMatrix(name='E', type='Full', nrow=2, ncol=2, labels=c('b0I', 'b0S', 'b1I', 'b1S'), values=0, free=TRUE)
xmat <- mxAlgebra( E %*% Z, name='x0')
modL <- list(amat, bmat, cmat, dmat, qmat, rmat, xmat,
	pmat, umat, zmat, emat)

# Create/Run model
indivmodels <- list()
for(k in 1:nSubj){
	dataL[[k]]$z <- z[k]
	DataSetForSubjectK <- dataL[[k]]
	indivmodels[[k]] <- mxModel(name=modNames[k],
		modL, expSS,
		mxFitFunctionML(),
		mxData(DataSetForSubjectK, type='raw')) 
}
multiSubjGrowthDef <- mxModel(name="MultiGrowthDef", indivmodels,
	mxFitFunctionMultigroup(modNames))

multiSubjGrowthDefRun <- mxRun(multiSubjGrowthDef)

# Inspect parameters
summary(multiSubjGrowthDefRun)

# Compare to linear mixed effects
require(lme4)
tallData <- do.call(rbind, dataL)
gd <- lmer(y ~ 1 + time*z + (1 + time | ID), data=tallData)

gdfix <- fixef(gd)
gdcov <- as.data.frame(VarCorr(gd))$vcov
cmpd <- cbind(StateSpace=coef(multiSubjGrowthDefRun), mixed=c(gdcov[c(4, 1, 3, 2)], gdfix))

omxCheckCloseEnough(cmpd[,1], cmpd[,2], 0.015)


#------------------------------------------------------------------------------
# Done
