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
# Date: 2014.03.13
# Filename: StateSpaceErrorCheck.R
# Purpose: Test the error checking for the state space expectation.
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# Revision History
# Thu 13 Mar 2014 12:57:28 Central Daylight Time -- Michael Hunter created file
# 


#--------------------------------------------------------------------
# Load required packages

require(OpenMx)
require(mvtnorm) # used to generate data



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
# Fit state space model to data via OpenMx package


smod <- mxModel(
	name='State Space Error Check',
	mxMatrix(name='A', values=tA, nrow=xdim, ncol=xdim),
	mxMatrix(name='B', values=0, nrow=xdim, ncol=udim, free=FALSE),
	mxMatrix(name='C', values=tC, nrow=ydim, ncol=xdim, free=(tC!=0), dimnames=list(rownames(ty), rownames(tx))),
	mxMatrix(name='D', values=0, nrow=ydim, ncol=udim, free=FALSE),
	mxMatrix(name='Q', type='Diag', values=diag(tQ), nrow=xdim, ncol=xdim, free=FALSE),
	mxMatrix(name='R', type='Diag', values=diag(tR), nrow=ydim, ncol=ydim, free=TRUE),
	mxMatrix(name='x', values=x0, nrow=xdim, ncol=1, free=FALSE),
	mxMatrix(name='P', values=P0, nrow=xdim, ncol=xdim, free=FALSE),
	mxMatrix("Zero", udim, 1, name="u"),
	mxData(observed=t(ty), type='raw'),
	mxExpectationStateSpace(A='A', B='B', C='C', D='D', Q='Q', R='R', x0='x', P0='P', u='u'),
	mxFitFunctionML()
)



serr <- mxModel(smod,
		mxMatrix(name='A', values=.8, nrow=xdim+1, ncol=xdim)
)
omxCheckError(mxRun(serr), "The A matrix is not the correct size in the state space expectation of model 'State Space Error Check'.  It is 4 by 3 and should be 3 by 3.")


serr <- mxModel(smod,
	mxMatrix(name='B', values=0, nrow=xdim+1, ncol=udim, free=FALSE)
)
omxCheckError(mxRun(serr), "The B matrix is not the correct size in the state space expectation of model 'State Space Error Check'.  It is 4 by 2 and should be 3 by 2.")



serr <- mxModel(smod,
	mxMatrix(name='B', values=0, nrow=xdim, ncol=udim+1, free=FALSE)
)
omxCheckError(mxRun(serr), "The D matrix is not the correct size in the state space expectation of model 'State Space Error Check'.  It is 9 by 2 and should be 9 by 3.")



serr <- mxModel(smod,
	mxMatrix(name='C', values=.4, nrow=ydim+1, ncol=xdim, dimnames=list(c(rownames(ty), "oops"), rownames(tx)))
)
omxCheckError(mxRun(serr), "The R matrix is not the correct size in the state space expectation of model 'State Space Error Check'.  It is 9 by 9 and should be 10 by 10.")



serr <- mxModel(smod,
	mxMatrix(name='C', values=.4, nrow=ydim, ncol=xdim+1, dimnames=list(rownames(ty), c(rownames(tx), "oops")))
)
omxCheckError(mxRun(serr), "The A matrix is not the correct size in the state space expectation of model 'State Space Error Check'.  It is 3 by 3 and should be 4 by 4.")



serr <- mxModel(smod,
	mxMatrix(name='D', values=0, nrow=ydim+1, ncol=udim, free=FALSE)
)
omxCheckError(mxRun(serr), "The D matrix is not the correct size in the state space expectation of model 'State Space Error Check'.  It is 10 by 2 and should be 9 by 2.")



serr <- mxModel(smod,
	mxMatrix(name='D', values=0, nrow=ydim, ncol=udim+1, free=FALSE)
)
omxCheckError(mxRun(serr), "The D matrix is not the correct size in the state space expectation of model 'State Space Error Check'.  It is 9 by 3 and should be 9 by 2.")



serr <- mxModel(smod,
	mxMatrix("Zero", udim+1, 1, name="u")
)
omxCheckError(mxRun(serr), "The u matrix is not the correct size in the state space expectation of model 'State Space Error Check'.  It is 3 by 1 and should be 2 by 1.")



serr <- mxModel(smod,
	mxMatrix("Zero", udim, 1+1, name="u")
)
omxCheckError(mxRun(serr), "The u matrix is not the correct size in the state space expectation of model 'State Space Error Check'.  It is 2 by 2 and should be 2 by 1.")


# Done

