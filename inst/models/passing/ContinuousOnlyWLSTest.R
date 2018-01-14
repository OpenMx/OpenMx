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
# Date: 2012.10.02
# Filename: wlsContTest.R
# Purpose: Test the mxDataWLS function for continuous-only data in an
#  OpenMx model.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Revision history
#  Wed 20 Aug 2014 13:25:04 Central Daylight Time -- Michael Hunter modified to use mxDataWLS.


#--------------------------------------
# Needed packages


require(OpenMx)
data(Bollen)

Bollen1 <- Bollen
Bollen1[1,'y1'] <- NA
omxCheckError(mxDataWLS(Bollen1[, 1:8]),
              "All continuous data with missingness cannot be handled in the WLS framework. Use na.omit(yourDataFrame) to remove rows with missing values or use maximum likelihood instead")

#--------------------------------------
# Set up  model matrices

manvar <- names(Bollen[, 1:8])

lval <- matrix(
		c(1, 0,
		  1, 0,
		  1, 0,
		  1, 0,
		  0, 1,
		  0, 1,
		  0, 1,
		  0, 1),
		byrow=TRUE,
		ncol=2, nrow=8)
lfre <- matrix(as.logical(lval), ncol=2)
lfre[1, 1] <- FALSE
lfre[5, 2] <- FALSE
llab <- matrix(c(paste("lam", 1:4, sep=""), rep(NA, 8), paste("lam", 1:4, sep="")), ncol=2)


lx <- mxMatrix(name="Lam", values=lval, free=lfre, ncol=2, nrow=8, labels=llab, dimnames=list(manvar, c("F1", "F2")))


td <- mxMatrix(name="Theta", type="Symm", ncol=8,
	values=
	c(.8,  0,  0,  0, .2,  0,  0,  0,
	      .8,  0, .2,  0, .2,  0,  0,
	          .8,  0,  0,  0, .2,  0,
	              .8,  0,  0,  0, .2,
	                  .8,  0,  0,  0,
	                      .8,  0, .2,
	                          .8,  0,
	                              .8),
	free=c(T,F,F,F,T,F,F,F,
	         T,F,T,F,T,F,F,
	           T,F,F,F,T,F,
	             T,F,F,F,T,
	               T,F,F,F,
	                 T,F,T,
	                   T,F,
	                     T),
	dimnames=list(manvar, manvar)
)
diag(td$labels) <- paste("var", 1:8, sep="")
selMat <- matrix(
	  c(5,1,
		4,2,
		6,2,
		7,3,
		8,4,
		8,6), ncol=2, byrow=TRUE)
td$labels[selMat] <- paste("cov", c(51, 42, 62, 73, 84, 86), sep="")
td$labels[selMat[,2:1]] <- paste("cov", c(51, 42, 62, 73, 84, 86), sep="")

ph <- mxMatrix(name="Phi", type="Symm", ncol=2, free=T, values=c(.8, .2, .8), labels=paste("phi", c(1, 12, 2), sep=""), dimnames=list(c("F1", "F2"), c("F1", "F2")))


#--------------------------------------
# Set-up WLS model

wlsMod <- mxModel("Test case for WLS Objective function from Bollen 1989",
	lx, ph, td,
	mxExpectationLISREL(LX=lx$name, PH=ph$name, TD=td$name),
	mxFitFunctionWLS(),
	mxDataWLS(Bollen[, 1:8])
)


# Run WLS model
wlsRun <- mxRun(wlsMod)
summary(wlsRun)
omxCheckTrue(is.null(wlsRun$output$calculatedHessian))

#TODO Fix summary for WLS data/fitfunctions
# Standard errors correct?
# observed statistics are not 0
# degrees of freedom should be correct when obs stats is fixed
# -2 log likelihood should be "Fit Function Value"?
# AIC, BIC for WLS?


# Compare parameter estimates
# WLS estimated parameters reported in Bollen (1989, p. 428, Table 9.4)
bollenParam <- c(l1=1.00, l2=1.11, l3=1.05, l4=1.16, e1=1.30, e2=7.12, e3=3.53, e4=2.72, e5=2.29, e6=3.77, e7=2.60, e8=3.45)

fitParam <- c(mxEval(Lam, wlsRun)[1:4,1], diag(mxEval(Theta, wlsRun)))

omxCheckCloseEnough(bollenParam, fitParam, epsilon=0.01)


#--------------------------------------
# Re-run with ML
mlMod <- mxModel(wlsMod, name="Rerun WLS model from Bollen as ML",
	mxFitFunctionML(),
	mxData(cov(Bollen[ , 1:8]), "cov", numObs=nrow(Bollen))
)
mlRun <- mxRun(mlMod)


#--------------------------------------
# Compare estimates

mlSum <- summary(mlRun)
wlsSum <- summary(wlsRun)

rms <- function(x, y){sqrt(mean((x-y)^2))}

# parameters are sort of close
omxCheckTrue(rms(mlSum$parameters[,5], wlsSum$parameters[,5]) < 0.7)

# standard errors are close
omxCheckTrue(rms(mlSum$parameters[,6], wlsSum$parameters[,6]) < 0.2)

# Chi square is on par
omxCheckWithinPercentError(mlSum$Chi, wlsSum$Chi, percent=85)

# Chi square df are the same
omxCheckEquals(mlSum$ChiDoF, wlsSum$ChiDoF)

# RMSEA confidence intervals overlap
omxCheckTrue( (wlsSum$RMSEA < mlSum$RMSEACI[2]) & (wlsSum$RMSEA > mlSum$RMSEACI[1]))
omxCheckTrue( (mlSum$RMSEA < wlsSum$RMSEACI[2]) & (mlSum$RMSEA >= wlsSum$RMSEACI[1]))

