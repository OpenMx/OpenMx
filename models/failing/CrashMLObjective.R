#
#   Copyright 2007-2010 The OpenMx Project
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

# -----------------------------------------------------------------------
# Program: MultivariateRegression_MatrixRaw.R  
#  Author: Ryne Estabrook
#    Date: 08 01 2009 
#
# Multivariate Regression model to estimate effect of independent on dependent variables
# Matrix style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
data(myRegDataRaw)

# add a column of data that is unused	
myRegDataRaw[,5] <- myRegDataRaw[,1]

covMatrix <- cov(myRegDataRaw[,c(1:4)])
#dimnames(covMatrix) <- list(colnames(myRegDataRaw)[1:4],colnames(myRegDataRaw)[1:4])
meansMatrix <- sapply(myRegDataRaw[,c(1:4)], mean)

#Create an MxModel object
# -----------------------------------------------------------------------
multivariateRegModel <- mxModel("Multiple Regression Matrix Specification", 
    mxData(
    	observed=covMatrix,
		means=meansMatrix,
    	type="cov",
		numObs=100
    ),
    mxMatrix(
    	type="Full", 
    	nrow=4, 
    	ncol=4,
        values=c(0,1,0,1,
                 0,0,0,0,
                 0,1,0,1,
                 0,0,0,0),
        free=c(F, T, F, T,
               F, F, F, F,
               F, T, F, T,
               F, F, F, F),
        labels=c(NA, "betawx", NA, "betawz",
                 NA,  NA,     NA,  NA, 
                 NA, "betayx", NA, "betayz",
                 NA,  NA,     NA,  NA),
        byrow=TRUE,
        name="A"
    ),
    mxMatrix(
    	type="Symm", 
    	nrow=4, 
    	ncol=4, 
        values=c(1,  0, 0,  0,
                 0,  1, 0, .5,
                 0,  0, 1,  0,
                 0, .5, 0,  1),
        free=c(T, F, F, F,
               F, T, F, T,
               F, F, T, F,
               F, T, F, T),
        labels=c("residualw",  NA,     NA,         NA,
                  NA,         "varx",  NA,        "covxz",
                  NA,          NA,    "residualy", NA,
                  NA,         "covxz", NA,        "varz"),
        byrow=TRUE,
        name="S"
    ),
    mxMatrix(
    	type="Iden",
    	nrow=4, 
    	ncol=4,
        name="F",
        dimnames=list(c("w", "x", "y", "z"),c("w", "x", "y", "z"))
    ),
    mxMatrix(
    	type="Full", 
    	nrow=1, 
    	ncol=4,
        values=c(0,0,0,0),
        free=c(T,T,T,T),
        labels=c("betaw","meanx","betay","meanz"),
        name="M"
    ),
	mxMatrix(type="Iden", nrow = 4, name = "I"),
	mxAlgebra(F %*% solve(I - A) %*% S %*% t(solve(I - A)) %*% t(F), name = 'covariance'),
	mxAlgebra(t(F %*% solve(I - A) %*% t(M)), name = 'means'),
    mxMLObjective("covariance", "means", dimnames=c("w", "x", "y", "z"))
)
      
multivariateRegFit<-mxRun(multivariateRegModel)

multivariateSummary <- summary(multivariateRegFit)


