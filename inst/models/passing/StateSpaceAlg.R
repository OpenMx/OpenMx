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
# Date: 2012.12.06
# Filename: StateSpaceOsc.R
# Purpose: Test the state space expectation.
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# Revision History
# Thu Dec 06 18:59:04 Central Standard Time 2012 -- Michael Hunter Checked in file to models/failing
# Thu 14 Feb 2013 15:52:57 Central Standard Time -- Michael Hunter realized the model actually worked.
# Thu 13 Feb 2014 15:06:58 Central Standard Time -- Michael Hunter removed mxConstraint and used parameter equal to result of mxAlgebra instead.
# 


#--------------------------------------------------------------------
# Load required packages

require(OpenMx)
require(mvtnorm) # used to generate data
#require(dlm) # only used if model is estimated with dlm for comparison


#--------------------------------------------------------------------
# Generate Data

xdim <- 1
udim <- 2
ydim <- 1
tdim <- 200
set.seed(950)
tA <- matrix(c(-.7), xdim, xdim)
tB <- matrix(c(0), xdim, udim)
tC <- matrix(c(1), ydim, xdim)
tD <- matrix(c(0), ydim, udim)
tQ <- matrix(c(0), xdim, xdim)
tR <- matrix(c(0), ydim, ydim); diag(tR) <- runif(ydim)

x0 <- matrix(c(rnorm(xdim)), xdim, 1)
P0 <- diag(c(runif(xdim)), nrow=xdim)
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

Astart <- tA

smod <- mxModel(
	name='StateSpaceExample',
	mxMatrix(name='A', values=Astart, nrow=xdim, ncol=xdim, free=TRUE, labels='a'),
	mxMatrix(name='B', values=0, nrow=xdim, ncol=udim, free=FALSE),
	mxMatrix(name='C', values=tC, nrow=ydim, ncol=xdim, free=FALSE, dimnames=list(rownames(ty), rownames(tx))),
	mxMatrix(name='D', values=0, nrow=ydim, ncol=udim, free=FALSE),
	# Note Factor error matrix is fixed!  This is for model identification.
	# I happen to fix the variances to their true values.
	mxMatrix(name='Q', type='Diag', values=diag(tQ), nrow=xdim, ncol=xdim, free=FALSE),
	mxMatrix(name='R', type='Diag', values=diag(tR), nrow=ydim, ncol=ydim, free=TRUE),
	mxMatrix(name='x', values=x0, nrow=xdim, ncol=1, free=FALSE),
	mxMatrix(name='P', values=P0, nrow=xdim, ncol=xdim, free=FALSE),
	mxMatrix("Zero", udim, 1, name="u"),
	mxData(observed=t(ty), type='raw'),
	mxExpectationStateSpace(A='A', B='B', C='C', D='D', Q='Q', R='R', x0='x', P0='P', u='u'),
	mxFitFunctionML(rowDiagnostics=TRUE)
)


# Uncomment for degugging
#smod <- mxOption(smod, 'Calculate Hessian', 'No')
#smod <- mxOption(smod, 'Standard Errors', 'No')
#smod <- mxOption(smod, 'Major iterations', 0)



srun <- mxRun(smod)

# Algebra works for A matrix
smodA <- mxModel(model=smod, name="State Space Example with Algebra for A Matrix",
	mxAlgebra(A, name="A2"),
	mxExpectationStateSpace(A='A2', B='B', C='C', D='D', Q='Q', R='R', x0='x', P0='P', u='u')
)

srunA <- mxRun(smodA)

omxCheckCloseEnough(summary(srun)$parameters[,5:6], summary(srunA)$parameters[,5:6], epsilon=1e-8)


# Algebra fails for C matrix
smodC <- mxModel(model=smod, name="State Space Example with Algebra for C Matrix",
	mxAlgebra(C, name="C2", dimnames=list(rownames(ty), rownames(tx))),
	mxExpectationStateSpace(A='A', B='B', C='C2', D='D', Q='Q', R='R', x0='x', P0='P', u='u')
)

srunC <- mxRun(smodC)

omxCheckCloseEnough(summary(srunC)$parameters[,5:6], summary(srunA)$parameters[,5:6], epsilon=1e-8)

# As a submodel
smodS <- mxModel(model="Sub", smod, mxFitFunctionAlgebra("StateSpaceExample.fitfunction"))
srunS <- mxRun(smodS)

omxCheckCloseEnough(summary(srun)$parameters[,5:6], summary(srunS)$parameters[,5:6], epsilon=1e-8)


# Check submodel likelihood evaluation
omxCheckCloseEnough(mxEval(Sub.fitfunction, srunS), 434.5657, epsilon=1e-3)

omxCheckCloseEnough(attr(mxEval(StateSpaceExample.fitfunction, srunS), "likelihoods"), attr(srunS$StateSpaceExample.fitfunction$result, "likelihoods"), epsilon=1e-8)





