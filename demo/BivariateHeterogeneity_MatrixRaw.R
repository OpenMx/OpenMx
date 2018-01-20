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


# -----------------------------------------------------------------------------
# Program: BivariateHeterogeneity_MatrixRaw.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Heterogeneity
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Bivariate Heterogeneity model to test differences in means and 
#      variances across multiple groups
#      Matrix style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------------

require(OpenMx)

require(MASS)
# Load Libraries
# -----------------------------------------------------------------------------

set.seed(200)
rs=.5
xy1 <- mvtnorm::rmvnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
xy1 <- xy1[, order(apply(xy1, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
# Note: Users do NOT have to re-order their data columns.  This is only to make data generation the same on different operating systems: to fix an inconsistency with the mvtnorm::rmvnorm function in the MASS package.
# group 1
# --------------------------------------
set.seed(200)
rs=.4
xy2 <- mvtnorm::rmvnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
xy2 <- xy2[, order(apply(xy2, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
# Note: Users do NOT have to re-order their data columns.  This is only to make data generation the same on different operating systems: to fix an inconsistency with the mvtnorm::rmvnorm function in the MASS package.
# group 2
# --------------------------------------
# Simulate Data
# -----------------------------------------------------------------------------


selVars <- c("X","Y")
summary(xy1)
cov(xy1)
dimnames(xy1) <- list(NULL, selVars)
summary(xy2)
cov(xy2)
dimnames(xy2) <- list(NULL, selVars)
# Print Descriptive Statistics
# -----------------------------------------------------------------------------

chol1        <- mxMatrix( type="Lower", nrow=2,  ncol=2, 
                          free=T, values=.5, labels=c("Ch11","Ch21","Ch31"), 
                          name="chol1" )
expCov1      <- mxAlgebra( expression=chol1 %*% t(chol1), name="expCov1" )
expMean1     <- mxMatrix( type="Full", nrow=1, ncol=2, 
                          free=T, values=c(0,0), labels=c("mX1","mY1"), 
                          name="expMean1" )
dataRaw1     <- mxData( xy1, type="raw" ) 
exp1         <- mxExpectationNormal( covariance="expCov1", means="expMean1", selVars)
funML        <- mxFitFunctionML()
model1       <- mxModel("group1", 
                         dataRaw1, chol1, expCov1, expMean1, exp1, funML)

chol2        <- mxMatrix( type="Lower", nrow=2,  ncol=2, 
                          free=T, values=.5, labels=c("Ch12","Ch22","Ch32"), 
                          name="chol2" )
expCov2      <- mxAlgebra( expression=chol2 %*% t(chol2), name="expCov2" )
expMean2     <- mxMatrix( type="Full", nrow=1, ncol=2, 
                          free=T, values=c(0,0), labels=c("mX2","mY2"), 
                          name="expMean2" )
dataRaw2     <- mxData( xy2, type="raw" ) 
exp2         <- mxExpectationNormal( covariance="expCov2", means="expMean2", selVars)
funML        <- mxFitFunctionML()
model2       <- mxModel("group2", 
                         dataRaw2, chol2, expCov2, expMean2, exp2, funML)

fun           <- mxFitFunctionMultigroup(c("group1.fitfunction", "group2.fitfunction"))
bivHetModel   <- mxModel("bivariate Heterogeneity Matrix Specification",
                        model1, model2, fun )

bivHetFit <- mxRun(bivHetModel)
expMean1Het <- mxEval(group1.expMean1, bivHetFit)
expMean2Het <- mxEval(group2.expMean2, bivHetFit)
expCov1Het <- mxEval(group1.expCov1, bivHetFit)
expCov2Het <- mxEval(group2.expCov2, bivHetFit)
LLHet <- bivHetFit$output$fit

expMean1Het; expMean2Het; expCov1Het; expCov2Het; LLHet

# Fit Heterogeneity Model
# -----------------------------------------------------------------------------

bivHomModel <- bivHetModel
bivHomModel[['group2.chol2']]$labels <- bivHomModel[['group1.chol1']]$labels
bivHomModel[['group2.expMean2']]$labels <- bivHomModel[['group1.expMean1']]$labels

bivHomFit <- mxRun(bivHomModel)
expMean1Hom <- mxEval(group1.expMean1, bivHomFit)
expMean2Hom <- mxEval(group2.expMean2, bivHomFit)
expCov1Hom <- mxEval(group1.expCov1, bivHomFit)
expCov2Hom <- mxEval(group2.expCov2, bivHomFit)
LLHom <- bivHomFit$output$fit

expMean1Hom; expMean2Hom; expCov1Hom; expCov2Hom; LLHom

Chi= LLHom-LLHet
LRT= rbind(LLHet,LLHom,Chi)
LRT
# Fit Homnogeneity Model
# -----------------------------------------------------------------------------

omxCheckCloseEnough(LLHet, 10927.4024, .001)
omxCheckCloseEnough(c(expCov1Het), c(1.0656, 0.4752, 0.4752, 0.9292), .001)
omxCheckCloseEnough(c(expMean1Het), c(0.0582, 0.0063), .001)
omxCheckCloseEnough(c(expCov2Het), c(1.0728, 0.3739, 0.3739, 0.9283),.001)
omxCheckCloseEnough(c(expMean2Het), c(0.0596, 0.0028),.001)

omxCheckCloseEnough(LLHom, 10936.95,.001)
omxCheckCloseEnough(c(expCov1Hom), c(1.0692, 0.4245, 0.4245, 0.9287),.001)
omxCheckCloseEnough(c(expMean1Hom), c(0.0589, 0.0046),.001)
omxCheckCloseEnough(c(expCov2Hom), c(1.0692, 0.4245, 0.4245, 0.9287),.001)
omxCheckCloseEnough(c(expMean2Hom), c(0.0589, 0.0046),.001)
# Compare OpenMx results to Mx results 
# (LL: likelihood; expCov: expected covariance, expMean: expected means)
# -----------------------------------------------------------------------------
