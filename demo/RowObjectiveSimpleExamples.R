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


#------------------------------------------------
# Author: Michael Hunter
# Filename: rowObjectiveSimpleExamples.R
# Purpose: Test the mxRowObjective Function
# Revision History
#  Mon Apr 11 20:54:30 EDT 2011 -- Created file
#  2011.06.17 -- Submitted file to SVN Repository
#------------------------------------------------


#------------------------------------------------
library(OpenMx)

#------------------------------------------------
# Model that adds two data columns row-wise, then sums that column
# Notice no optimization is performed here.
set.seed(159)
xdat <- data.frame(a=rnorm(10, mean=4.2), b=1:10) # Make data set

amod <- mxModel(	
	name='Row Model',
	mxData(observed=xdat, type='raw'),
	mxAlgebra(sum(filteredDataRow), name = 'rowAlgebra'),
	mxAlgebra(sum(rowResults), name = 'reduceAlgebra'),
	mxRowObjective(
		rowAlgebra='rowAlgebra',
		reduceAlgebra='reduceAlgebra',
		dimnames=c('a','b'))
)

amodFit <- mxRun(amod)
mxEval(objective, amodFit)
sum(xdat)

omxCheckCloseEnough(mxEval(objective, amodFit), sum(xdat), epsilon=10^(-5))


#------------------------------------------------
# Model that finds the parameters that minimizes the sum of the
#  squared difference between the parameters and a data rows.
# This is a least squares estimation of the means of the columns.
xdat

bmod <- mxModel(
	name='Estimation Row Model',
	mxData(observed=xdat, type='raw'),
	mxMatrix(values=.75, ncol=2, nrow=1, free=TRUE, name='M'),
	mxAlgebra((filteredDataRow-M)%^%2, name='rowAlgebra'),
	mxAlgebra(sum(rowResults), name='reduceAlgebra'),
	mxRowObjective(
		rowAlgebra='rowAlgebra',
		reduceAlgebra='reduceAlgebra',
		dimnames=c('a', 'b'))
)

bmodFit <- mxRun(bmod)
bmodFit$M@values
colMeans(xdat)

omxCheckCloseEnough(as.vector(mxEval(M, bmodFit)), as.vector(colMeans(xdat)), epsilon=10^(-5))



#------------------------------------------------
# Model that finds the parameters that minimizes the sum of the
#  squared difference between the parameters and a data rows.
# This is a least squares estimation of the means of the columns,
#  taking into account missingness.
# This script fails even when the NAs are removed.
#  It seems to be a problem with the omxSelect*.

# Comment out the next two lines to test the use of 
#  omxSelectCols when there is no missingness.
xdat$a[3] <- NA
xdat$b[5] <- NA
xdat


cmod <- mxModel(
	name='Estimation Row Model with Missingness',
	mxData(observed=xdat, type='raw'),
	mxMatrix(values=.75, ncol=2, nrow=1, free=TRUE, name='M'),
	mxAlgebra(omxSelectCols(M, existenceVector), name='fM'),
	mxAlgebra((filteredDataRow-fM)%^%2, name='rowAlgebra'),
	mxAlgebra(sum(rowResults), name='reduceAlgebra'),
	mxRowObjective(
		rowAlgebra='rowAlgebra',
		reduceAlgebra='reduceAlgebra',
		dimnames=c('a', 'b'))
)

cmodFit <- mxRun(cmod)
cmodFit$M@values
colMeans(xdat, na.rm=T)

omxCheckCloseEnough(as.vector(mxEval(M, cmodFit)), as.vector(colMeans(xdat, na.rm=T)), epsilon=10^(-5))




#------------------------------------------------
set.seed(135)
nobs <- 13
adat <- data.frame(x=rnorm(nobs))

dmod <- mxModel(
	name='I will run fast on OpenMx',
	mxMatrix(name='A', nrow=nobs, ncol=1, free=T, values=0.1),
	mxMatrix(name='X', nrow=nobs, ncol=1, free=F, values=as.matrix(adat)),
	mxAlgebra((X-A) %^% 2, name='Row'),
	mxAlgebra(sum(Row), name='Red'),
	mxAlgebraObjective('Red')
)

dmodRun <- mxRun(dmod) # runs super fast := 0.07 sec
omxCheckCloseEnough(mxEval(A, dmodRun), as.matrix(adat), epsilon=10^(-5))



#------------------------------------------------
robj1 <- function(model, state) {
	a <- model$A@values
	x <- model$X@values
	return(sum((x - a) ^ 2))
}

emod <- mxModel(
	name='I will run slow on OpenMx',
	mxMatrix(name='A', nrow=nobs, ncol=1, free=T, values=0.1),
	mxMatrix(name='X', nrow=nobs, ncol=1, free=F, values=as.matrix(adat)),
	mxRObjective(robj1)
)

emodRun <- mxRun(emod) # runs super slow := 10.5 sec
omxCheckCloseEnough(mxEval(A, emodRun), as.matrix(adat), epsilon=10^(-5))


