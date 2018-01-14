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

library(mvtnorm)

#--------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2012.12.06
# Filename: StateSpaceOsc.R
# Purpose: Test the state space expectation.
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# Revision History
# Thu Dec 06 18:59:04 Central Standard Time 2012 -- Michael Hunter Checked in file to models/failing
# Thu 14 Feb 2013 15:52:57 Central Standard Time -- Michael Hunter realized the model actually worked.
# 


#--------------------------------------------------------------------
# Load required packages

require(OpenMx)
require(mvtnorm) # used to generate data
#require(dlm) # only used if model is estimated with dlm for comparison


#--------------------------------------------------------------------
# Generate Data

xdim <- 3
udim <- 2
ydim <- 9
tdim <- 200
set.seed(948)
tA <- matrix(c(-.4, 0, 0, 0, -.9, .1, 0, -.1, -.9), xdim, xdim)
tB <- matrix(c(0), xdim, udim)
tC <- matrix(c(runif(3, .4, 1), rep(0, ydim), runif(3, .4, 1), rep(0, ydim), runif(3, .4, 1)), ydim, xdim)
tD <- matrix(c(0), ydim, udim)
tQ <- matrix(c(0), xdim, xdim); diag(tQ) <- runif(xdim)
tR <- matrix(c(0), ydim, ydim); diag(tR) <- runif(ydim)

x0 <- matrix(c(rnorm(xdim)), xdim, 1)
P0 <- diag(c(runif(xdim)))
tx <- matrix(0, xdim, tdim+1)
tu <- matrix(0, udim, tdim)
ty <- matrix(0, ydim, tdim)

tx[,1] <- x0
for(i in 2:(tdim+1)){
	tx[,i] <- tA %*% tx[,i-1] + tB %*% tu[,i-1] + t(rmvnorm(1, rep(0, xdim), tQ))
	ty[,i-1] <- tC %*% tx[,i-1] + tD %*% tu[,i-1] + t(rmvnorm(1, rep(0, ydim), tR))
}

#plot(tx[1,], type='l')

rownames(ty) <- paste('y', 1:ydim, sep='')
rownames(tx) <- paste('x', 1:xdim, sep='')


#--------------------------------------------------------------------
# Fit state space model to data via dlm package
# For posterity show how the same model would be estimated in the dlm package.
# This is how the values I validated the estimation for OpenMx,
#  i.e. by comparing the estimates from dlm and OpenMx.
# Note that in my (mhunter) experience OpenMx is much faster (25x in this example).

#mfun <- function(x){
#	mG <- matrix(c(x[1], 0, 0, 0, x[2], x[3], 0, -x[3], x[2]), xdim, xdim)
#	mW <- tQ # diag(x[4:6])
#	mF <- matrix(c(x[7:9], rep(0, ydim), x[10:12], rep(0, ydim), x[13:15]), ydim, xdim)
#	mV <- diag(x[16:24])
#	mM <- x0
#	mC <- P0
#	return(dlm(FF=mF, V=mV, GG=mG, W=mW, m0=mM, C0=mC))
#}


#tinit <- c(-.4, -.9, .1, diag(tQ), tC[tC!=0], diag(tR))
#mfun(tinit)

#dlmBegin <- Sys.time()
#mfit <- dlmMLE(y=t(ty), parm=tinit, build=mfun, lower=c(rep(NA, 3), rep(0.00001, 3), rep(NA, 9), rep(0.00001, 9)), control=list(maxit=200))
#dlmEnd <- Sys.time()
#mfun(mfit$par)

#mfun(mfit$par)$GG
#tA

#mfun(mfit$par)$FF
#tC

#diag(mfun(mfit$par)$W)
#diag(tQ)

#diag(mfun(mfit$par)$V)
#diag(tR)



#--------------------------------------------------------------------
# Fit state space model to data via OpenMx package

smod <- mxModel(
	name='StateSpaceExample',
	mxMatrix(name='A', values=tA, nrow=xdim, ncol=xdim, free=c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), labels=c('a', NA, NA, NA, 'b', 'c', NA, 'd', 'b'), ubound=c(-.6, rep(NA, 8))),
	mxAlgebra(name='csym', -c),
	mxConstraint(name='ccon', d == csym),
	mxMatrix(name='B', values=0, nrow=xdim, ncol=udim, free=FALSE),
	mxMatrix(name='C', values=tC, nrow=ydim, ncol=xdim, free=(tC!=0), dimnames=list(rownames(ty), rownames(tx))),
	mxMatrix(name='D', values=0, nrow=ydim, ncol=udim, free=FALSE),
	# Note Factor error matrix is fixed!  This is for model identification.
	# I happen to fix the variances to their true values.
	mxMatrix(name='Q', type='Diag', values=diag(tQ), nrow=xdim, ncol=xdim, free=FALSE),
	mxMatrix(name='R', type='Diag', values=diag(tR), nrow=ydim, ncol=ydim, free=TRUE),
	mxMatrix(name='x', values=x0, nrow=xdim, ncol=1, free=FALSE),
	mxMatrix(name='P', values=P0, nrow=xdim, ncol=xdim, free=FALSE),
	mxMatrix("Zero", udim, 1, name="u"),
	mxExpectationStateSpace(A='A', B='B', C='C', D='D', Q='Q', R='R', x='x', P='P', u='u'),
	mxFitFunctionML()
)

container <- mxModel("Container", smod,
                     mxData(observed=t(ty)[1:200,], type='raw'),
                     mxFitFunctionMultigroup(groups="StateSpaceExample.fitfunction"))

# Uncomment for degugging
#smod <- mxOption(smod, 'Calculate Hessian', 'No')
#smod <- mxOption(smod, 'Standard Errors', 'No')
#smod <- mxOption(smod, 'Major iterations', 0)


#ssmBegin <- Sys.time()
container <- mxRun(container)
#ssmEnd <- Sys.time()

#ssmEnd-ssmBegin
#dlmEnd-dlmBegin
# OpenMx is 24.6 times faster then dlm

# Check likelihoods of initial parameters
# -2LL
#summary(srun)$Minus2LogLikelihood # when major iterations is 0
#2*dlmLL(y=t(ty), mod=mfun(tinit)) + 200*9*log(2*pi) # dlm gives back -LL - CONST, so adjust it.  200 is N, 9 is k

srun <- container$submodels$"StateSpaceExample"$matrices

dlmEstA <- matrix(c(
	-0.7911864,  0.0000000,  0.0000000,
	 0.0000000, -0.8960419, -0.1064521,
	 0.0000000,  0.1064521, -0.8960419),
	3, 3, byrow=TRUE)

dlmEstC <- c( #nonzero factor loadings
	0.1798166, 0.4718692, 0.4547457,
	0.9500226, 0.6864060, 0.9554287,
	0.4287551, 0.9650694, 0.4956449)

dlmEstR <- c( #diagonal manifest error cov
	0.3798346, 0.8068893, 1.0383961,
	0.6261729, 0.1331556, 0.7761499,
	0.7538665, 0.9816791, 1.2166798)



omxCheckCloseEnough(srun$A$values, dlmEstA, epsilon=0.001)
omxCheckCloseEnough(srun$C$values[srun$C$free], dlmEstC, epsilon=0.001)
omxCheckCloseEnough(diag(srun$R$values), dlmEstR, epsilon=0.001)


