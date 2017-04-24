#
#   Copyright 2007-2017 The OpenMx Project
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

#   Script (by Rob K.) that demonstrates use of the bivariate-lognormal likelihood function, via 
#   mxFitFunctionRow().  The likelihood is actually that of a bivariate lognormal distribution shifted 
#   downwardly by 1.  Fitting this distribution is equivalent to applying a log(x+1) transformation to both
#   variables, and then applying the bivariate-normal likelihood to the data.  The "plus one" part is so that
#   data values of zero (the smallest they can be in this example) are made positive before trying to take 
#   their natural logarithm.

#   Ordinarily, it is not valid to compare the AIC of a model fitted to transformed data to that of a model 
#   fitted to the untransformed data.  However, incorporating the transformation into the likelihood as done
#   here would allow for valid comparison of this model's AIC to that of another model which assumes a different
#   continuous distribution (e.g., the default bivariate normal).

#   This script is being used as a regression test, to check for runtime problems with row fitfunctions that 
#   result from missing data.


library(OpenMx)
set.seed(1234)
x <- rpois(500,1)
y <- rpois(500,1)
y[500] <- NA
varNames <- c('x','y')
dat <- cbind(x,y,weight=1.0)

mymod <- mxModel(
	"bivlognorm",
	mxData(dat, "raw", weight = "weight"),
	mxMatrix(type="Full",nrow=1,ncol=2,free=T,values=0.1,labels="m",name="Mu",
					 dimnames=list(NULL,varNames)),
	mxMatrix(type="Symm",nrow=2,ncol=2,free=T,values=c(0.4,0,0.4),labels=c("v","c","v"),name="Sigma",
					 lbound=c(0.0001,-Inf,0.0001),
					 dimnames=list(varNames,varNames)),
	mxMatrix(type="Unit",nrow=1,ncol=2,name="ONE"),
	mxAlgebra(
		2*log(2*3.1415927*prod(filteredDataRow + omxSelectCols(ONE,existenceVector))) +
			log(det(omxSelectRowsAndCols(Sigma,existenceVector))) + 
			( (log(filteredDataRow+omxSelectCols(ONE,existenceVector))-omxSelectCols(Mu,existenceVector)) %*% 
					solve(omxSelectRowsAndCols(Sigma,existenceVector)) %*% 
					t(log(filteredDataRow+omxSelectCols(ONE,existenceVector))-omxSelectCols(Mu,existenceVector)) ),
		name="unweightedRowAlgebra"),
	mxAlgebra(data.weight * unweightedRowAlgebra, name="rowAlgebra"),
	mxAlgebra(sum(rowResults), name="reduceAlgebra"),
	mxFitFunctionRow(rowAlgebra='rowAlgebra',
								 reduceAlgebra='reduceAlgebra',
								 dimnames=c("x","y"))
)
myrun <- mxRun(mymod)
omxCheckCloseEnough(myrun$output$fit, 2644.072, .01)

boot <- mxBootstrap(myrun, 10)
bq1 <- summary(boot)[["bootstrapQuantile"]]
omxCheckTrue(all(apply(bq1,1,diff) > 0))
