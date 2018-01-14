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


#--------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2014.02.17
# Filename: StateSpaceInputs.R
# Purpose: Test the state space expectation with Inputs.
#  Note that the estimates here are only compared to the previous
#  estimates produced by this same program.  So the ground truth is
#  somewhat circular.  However, the estimates correspond quite nicely
#  to the true generating parameters.  In fact, these estimates are
#  even better than those of the same model without inputs.
#  Cf. StateSpaceOsc.R
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# Revision History
# Thu Dec 06 18:59:04 Central Standard Time 2012 -- Michael Hunter Checked in file to models/failing
# Thu 14 Feb 2013 15:52:57 Central Standard Time -- Michael Hunter realized the model actually worked.
# Thu 13 Feb 2014 15:06:58 Central Standard Time -- Michael Hunter removed mxConstraint and used parameter equal to result of mxAlgebra instead.
# Mon 17 Feb 2014 19:30:07 Central Standard Time -- Michael Hunter added inputs and created file from StateSpaceOsc.R
# Thu 20 Mar 2014 16:00:34 Central Daylight Time -- Michael Hunter added estimated value checking against previous estimates.
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
set.seed(2227)
tA <- matrix(c(-.4, 0, 0, 0, -.9, .1, 0, -.1, -.9), xdim, xdim)
tB <- matrix(c(3.1, 2.7, -2.5, 0, 0, 0), xdim, udim)
tC <- matrix(c(runif(3, .4, 1), rep(0, ydim), runif(3, .4, 1), rep(0, ydim), runif(3, .4, 1)), ydim, xdim)
tD <- matrix(c(rep(0, ydim), rep(.2, ydim)), ydim, udim)
tQ <- matrix(c(0), xdim, xdim); diag(tQ) <- runif(xdim)
tR <- matrix(c(0), ydim, ydim); diag(tR) <- runif(ydim)

x0 <- matrix(c(rnorm(xdim)), xdim, 1)
P0 <- diag(c(runif(xdim)))
tx <- matrix(0, xdim, tdim+1)
#tu <- matrix(rnorm(udim*tdim), udim, tdim) # Note: Random inputs!
tu <- matrix(c(rep(1, tdim), rnorm(tdim)), udim, tdim, byrow=TRUE) # Note: Constant and random inputs
ty <- matrix(0, ydim, tdim)

tx[,1] <- x0
#for(i in 1:tdim){
#	tx[,i+1] <- tA %*% tx[,i] + tB %*% tu[,i] + t(rmvnorm(1, rep(0, xdim), tQ))
#	ty[,i] <- tC %*% tx[,i] + tD %*% tu[,i] + t(rmvnorm(1, rep(0, ydim), tR))
#}
for(i in 2:(tdim+1)){
	tx[,i] <- tA %*% tx[,i-1] + tB %*% tu[,i-1] + t(rmvnorm(1, rep(0, xdim), tQ))
	ty[,i-1] <- tC %*% tx[,i] + tD %*% tu[,i-1] + t(rmvnorm(1, rep(0, ydim), tR))
}

#plot(tx[1,], type='l')

rownames(ty) <- paste('y', 1:ydim, sep='')
rownames(tu) <- paste('u', 1:udim, sep='')
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

AStart <- matrix(0, xdim, xdim)
AStart[tA!=0] <- .5
BStart <- c(rep(1, xdim), rep(0, xdim))
CStart <- matrix(0, ydim, xdim)
CStart[tC!=0] <- .8
DStart <- c(rep(0, ydim), rep(.5, ydim))
RStart <- .2


smod <- mxModel(
	name='State Space Example with Inputs',
	mxMatrix(name='A', values=AStart, nrow=xdim, ncol=xdim, free=c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE), labels=c('a', NA, NA, NA, 'b', 'c', NA, 'csym[1,1]', 'b')),
	mxAlgebra(name='csym', -c),
	mxMatrix(name='B', values=BStart, nrow=xdim, ncol=udim, free=c(rep(TRUE, xdim), rep(FALSE, xdim))),
	mxMatrix(name='C', values=CStart, nrow=ydim, ncol=xdim, free=(tC!=0), lbound=1e-6, dimnames=list(rownames(ty), rownames(tx))),
	mxMatrix(name='D', values=DStart, nrow=ydim, ncol=udim, free=c(rep(FALSE, ydim), rep(TRUE, ydim))),
	# Note Factor error matrix is fixed!  This is for model identification.
	# I happen to fix the variances to their true values.
	mxMatrix(name='Q', type='Diag', values=diag(tQ), nrow=xdim, ncol=xdim, free=FALSE),
	mxMatrix(name='R', type='Diag', values=RStart, nrow=ydim, ncol=ydim, free=TRUE),
	mxMatrix(name='x', values=x0, nrow=xdim, ncol=1, free=FALSE),
	mxMatrix(name='P', values=P0, nrow=xdim, ncol=xdim, free=FALSE),
	mxMatrix(type="Full", udim, 1, labels=c("data.u1", "data.u2"), name="u"),
	mxData(observed=cbind(t(ty), t(tu)), type='raw'),
	mxExpectationStateSpace(A='A', B='B', C='C', D='D', Q='Q', R='R', x0='x', P0='P', u='u', scores=TRUE),
	mxFitFunctionML()
)




srun <- mxRun(smod)


# Notice that the estimated parameters are close their true generating values
srun$A$values
tA


srun$B$values
tB

srun$C$values
tC


srun$D$values
tD



summary(srun)

prevEstA <- matrix(c(
	 -0.5120488,  0.0000000,  0.0000000,
	  0.0000000, -0.9138461, -0.0804081,
	  0.0000000,  0.0804081, -0.9138461),
	3, 3, byrow=TRUE)

prevEstB <- matrix(c(
	  3.586831,    0,
	  2.566541,    0,
	 -2.467185,    0),
	3, 2, byrow=TRUE
)

prevEstC <- c( #nonzero factor loadings
	0.8098549, 0.7388864, 0.8933953,
	0.6473688, 0.9509194, 0.7705916,
	0.9083092, 0.4261629, 0.8682982)

prevEstD <- c( #nonzero part of feedthrough matrix
	0.17368033, 0.23882681, 0.23983565,
	0.23048017, 0.26401039, 0.24613672,
	0.24795466, 0.17059682, 0.09830594)

prevEstR <- c( #diagonal manifest error cov
	0.45646303, 0.16792211, 0.54165256,
	0.07453542, 0.41766447, 0.33869364,
	0.21551307, 0.69513482, 0.69142726)



omxCheckCloseEnough(srun$A$values, prevEstA, epsilon=0.001)
omxCheckCloseEnough(srun$B$values, prevEstB, epsilon=0.01)
omxCheckCloseEnough(srun$C$values[srun$C$free], prevEstC, epsilon=0.001)
omxCheckCloseEnough(srun$D$values[srun$D$free], prevEstD, epsilon=0.001)
omxCheckCloseEnough(diag(srun$R$values), prevEstR, epsilon=0.001)

# Check to make sure definition variables can be evaluated
uAtRow10 <- mxEval(u, smod, compute=TRUE, defvar.row=10)
dAtRow10 <- mxEval(rbind(data.u1, data.u2), smod, compute=TRUE, defvar.row=10)
urAtRow10 <- mxEval(u, srun, compute=TRUE, defvar.row=10)
drAtRow10 <- mxEval(rbind(data.u1, data.u2), srun, compute=TRUE, defvar.row=10)

omxCheckCloseEnough(uAtRow10, dAtRow10, epsilon=1e-10)
omxCheckCloseEnough(uAtRow10, drAtRow10, epsilon=1e-10)
omxCheckCloseEnough(urAtRow10, drAtRow10, epsilon=1e-10)



# Check that OpenMx errors appropriately for a missing definition variable.
naData <- cbind(t(ty), t(tu))
naData[10,11] <- NA
naModel <- mxModel(smod, name='missingDef', mxData(observed=naData, type='raw'))

omxCheckError(naRun <- mxRun(naModel), "missingDef.data: NA in definition variable 'u2' row 10")

# Check for proper error message when there exist missing definition variables.


cor(t(tx), srun$expectation$xPredicted)
cor(t(tx), srun$expectation$xUpdated)
cor(t(tx), srun$expectation$xSmoothed)



