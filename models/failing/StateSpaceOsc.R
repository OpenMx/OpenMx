#
#   Copyright 2007-2012 The OpenMx Project
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
# Date: 2012.12.06
# Filename: StateSpaceOsc.R
# Purpose: Test the state space expectation.
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# Revision History
# Thu Dec 06 18:59:04 Central Standard Time 2012 -- Michael Hunter Checked in file to models/failing
# 


#--------------------------------------------------------------------


require(OpenMx)
require(mvtnorm)
#require(dlm)


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

rownames(ty) <- paste('y', 1:ydim)
rownames(tx) <- paste('x', 1:xdim)


#--------------------------------------------------------------------
# Fit state space model to data via dlm package

#mfun <- function(x){
#	mG <- matrix(c(x[1], 0, 0, 0, x[2], x[3], 0, -x[3], x[2]), xdim, xdim)
#	mW <- diag(x[4:6])
#	mF <- matrix(c(x[7:9], rep(0, ydim), x[10:12], rep(0, ydim), x[13:15]), ydim, xdim)
#	mV <- diag(x[16:24])
#	mM <- x0
#	mC <- P0
#	return(dlm(FF=mF, V=mV, GG=mG, W=mW, m0=mM, C0=mC))
#}

#tinit <- c(-.5, .9, .01, diag(tQ), tC[tC!=0], diag(tR))
#mfun(tinit)

#mfit <- dlmMLE(y=t(ty), parm=tinit, build=mfun, lower=c(rep(NA, 3), rep(0.00001, 3), rep(NA, 9), rep(0.00001, 9)), control=list(maxit=20))
#mfun(mfit$par)

#mfun(mfit$par)$GG
tA

#mfun(mfit$par)$FF
tC

#diag(mfun(mfit$par)$W)
diag(tQ)

#diag(mfun(mfit$par)$V)
diag(tR)

(Sinv <- solve(tC %*% (tA %*% P0 %*% t(tA) + tQ) %*% t(tC) + tR)) # S^-1

(Y <- tC %*% (tA %*% P0 %*% t(tA) + tQ)) # Y

(K <- t(Y) %*% Sinv) # Correct K

(Kin <- matrix(Y, nrow=3) %*% Sinv) # Incorrect K

# P = P - K C P
(tA %*% P0 %*% t(tA) + tQ) - K %*% tC %*% (tA %*% P0 %*% t(tA) + tQ)

(tA %*% P0 %*% t(tA) + tQ) - Kin %*% tC %*% (tA %*% P0 %*% t(tA) + tQ)

#--------------------------------------------------------------------
# Fit state space model to data via OpenMx package

smod <- mxModel(
	name='State Space Example',
	mxMatrix(name='A', values=tA, nrow=xdim, ncol=xdim, free=c(T, F, F, F, T, T, F, F, T), labels=c('a', NA, NA, NA, 'b', 'c', NA, 'd', 'b')),
	#mxAlgebra(name='csym', -c),
	#mxConstraint(name='ccon', d == csym),
	mxMatrix(name='B', values=0, nrow=xdim, ncol=udim, free=FALSE),
	mxMatrix(name='C', values=tC, nrow=ydim, ncol=xdim, free=(tC!=0), dimnames=list(rownames(ty), rownames(tx))),
	mxMatrix(name='D', values=0, nrow=ydim, ncol=udim, free=FALSE),
	mxMatrix(name='Q', type='Diag', values=diag(tQ), nrow=xdim, ncol=xdim, free=TRUE),
	mxMatrix(name='R', type='Diag', values=diag(tR), nrow=ydim, ncol=ydim, free=TRUE),
	mxMatrix(name='x', values=x0, nrow=xdim, ncol=1, free=FALSE),
	mxMatrix(name='P', values=P0, nrow=xdim, ncol=xdim, free=FALSE),
	mxData(observed=t(ty), type='raw'),
	imxExpectationStateSpace(A='A', B='B', C='C', D='D', Q='Q', R='R', x='x', P='P'),
	mxFitFunctionML()
)

smod <- mxOption(smod, 'Calculate Hessian', 'No')
smod <- mxOption(smod, 'Standard Errors', 'No')
smod <- mxOption(smod, 'Major iterations', 2)
srun <- mxRun(smod, onlyFrontend=FALSE)




