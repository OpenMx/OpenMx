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


#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2014.04.12
# Filename: MultilevelStateSpaceEx5.R
# Purpose: Replicate Mplus results reported by Gu et al. 2014 for a 2-level
#  structural equation model.  It's a factor model with "random intercepts".
#------------------------------------------------------------------------------


# Revision History
#  Mon 21 Apr 2014 07:23:58 Central Daylight Time -- Michael Hunter added comments for file distribution to
#  Steve Boker and Colleagues for Multilevel BG collaboration
#  


#------------------------------------------------------------------------------
# Load Package, read in data, set a variable or two

require(OpenMx)  #Note: must be r3345 or later

# Source of data
#http://www.infoagepub.com/products/content/files/serlinfiles_multilevel/Chapter%2010/Mplus/

# setwd("../Desktop")
# setwd("C:/Users/mhunter/Documents/Projects/OpenMx/MultilevelSSM")
# requires internet connection
data <- read.table("http://www.infoagepub.com/products/content/files/serlinfiles_multilevel/Chapter%2010/Mplus/Ch10Mplus.dat")

dataL <- list()
for(k in unique(data$V1)){
	dataL[[k]] <- data[data$V1==k,]
}

nx <- 3
ny <- 3
nu <- 1

#mxOption(NULL, "Default optimizer", "NPSOL")


#------------------------------------------------------------------------------
# Specify and fit the model


# The WITHIN model
LW <- 	mxMatrix(name="LamW", "Full", ny, nu, values=1, labels=paste('lw', 1:ny, sep=''), free=c(FALSE, TRUE, TRUE))
PW <- 	mxMatrix(name="PsiW", "Full", nu, nu, values=1, free=TRUE, labels="psw")
TW <- 	mxMatrix(name="ThdW", "Diag", ny, ny, values=1, free=TRUE, labels=paste('thw', 1:ny, sep=''))
CW <- 	mxAlgebra(name="CovW", LamW %*% PsiW %*% t(LamW) + ThdW)


# The BETWEEN model
LB <- 	mxMatrix(name="LamB", "Full", nx, nu, values=1, labels=paste('lb', 1:ny, sep=''), free=c(FALSE, TRUE, TRUE))
PB <- 	mxMatrix(name="PsiB", "Full", nu, nu, values=1, free=TRUE, labels="psb")
TB <- 	mxMatrix(name="ThdB", "Diag", nx, nx, values=1, free=TRUE, labels=paste('thb', 1:ny, sep=''))
CB <- 	mxAlgebra(name="CovB", LamB %*% PsiB %*% t(LamB) + ThdB)


# The state space matrices
a <- mxMatrix(name="A", "Diag", nx, nx, values=1)
b <- mxMatrix(name="B", "Zero", nx, nu)
c <- mxMatrix(name="C", "Diag", ny, nx, values=1, dimnames=list(c('V2', 'V3', 'V4'), c('u_1', 'u_2', 'u_3')))
d <- mxMatrix(name="D", "Full", ny, nu, values=1, free=TRUE, labels=c('g100', 'g200', 'g300'))
q <- mxMatrix(name="Q", "Zero", nx, nx)
u <- mxMatrix(name="u", "Full", nu, 1, values=1)
x0 <- mxMatrix(name="x0", "Zero", nx, 1)

# Create a list of mxModels, one for each level 2 unit (classroom)
indivmodels <- list()
for(k in unique(data$V1)){
	indivmodels[[k]] <- mxModel(name=paste('indiv', k, sep=''),
		a, b, c, d, q, x0, u,
		LW, PW, TW, CW, LB, PB, TB, CB,
		mxExpectationStateSpace(A="A", B="B", C="C", D="D", Q="Q", R="CovW", x0="x0", P0="CovB", u="u"),
		# note the R and P0 matrices
		mxFitFunctionML(),
		mxData(dataL[[k]], type='raw')) # note separate data
}

# There should be a slicker way to do this, but I haven't found it yet.
# I'd like to have one model where the SSM matrices live, one model for the WITHIN matrices,
#  one model for the BETWEEN matrices, and then a list of models that refer to these, each with
#  their own data.  I guess that's what flattening is doing.
# Still I'd like to use built-in expectation functions to create the within/between
#  expected covariance matrices.

groupfits <- paste("indiv", unique(data$V1), ".fitfunction", sep="")

m3s <- mxModel(model='Container', indivmodels, mxFitFunctionMultigroup(groupfits))
m3sRun <- mxRun(m3s)

summary(m3sRun) #runs in approx 60 seconds on my machine

# matches table 5
mplusParam <- c(3.677, 3.731, 3.706, .978, 1.070, 1.248, 1.769, 1.868, 1.510, 1.039, 1.091, .107, .025, .058, .037)
omxCheckCloseEnough(omxGetParameters(m3sRun), mplusParam, 0.001)

#Also check standard errors
mplusSE <- c(.073, .082, .080, .040, .044, .074, .063, .063, .065, .234, .228, .042, .019, .026, .023)
rd <- (summary(m3sRun)$parameters[,6] - mplusSE) / mplusSE
omxCheckCloseEnough(rd, rep(0, length(mplusSE)), 0.01)


#------------------------------------------------------------------------------
# Different way to do the same multilevel model.
# This way avoids the multigroup overhead.
# Instead of using multiple groups, it just re-initializes

# This example only re-initialized the covariances, bc the
# between cluster means are always zero.


# Create continue cluster variable
# 0 at the beginning of each cluster, 1 everywhere else.
dataC <- data
dataC$cclus <- c(0, 1-diff(dataC$V1))


a2 <- mxMatrix(name='A', 'Diag', nx, nx, labels='data.cclus')
sw <- mxMatrix(name='con', 'Full', 1, 1, labels='data.cclus')
q2 <- mxAlgebra((1-con) %x% CovB, name='Q')

simpleModel <- mxModel('SimpleMultilevel',
	a2, b, c, d, q2, x0, u, sw,
	LW, PW, TW, CW, LB, PB, TB, CB,
	mxExpectationStateSpace(A="A", B="B", C="C", D="D", Q="Q", R="CovW", x0="x0", P0="CovB", u="u"),
	# note the R and P0 matrices
	mxFitFunctionML(),
	mxData(dataC, type='raw')
)

simpleRun <- mxRun(simpleModel)

omxCheckCloseEnough(omxGetParameters(simpleRun), mplusParam, 0.001)

rd <- (summary(simpleRun)$parameters[,6] - mplusSE) / mplusSE
omxCheckCloseEnough(rd, rep(0, length(mplusSE)), 0.05)



