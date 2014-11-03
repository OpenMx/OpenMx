#
#   Copyright 2007-2014 The OpenMx Project
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
xy1 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
xy1 <- xy1[, order(apply(xy1, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
# Note: Users do NOT have to re-order their data columns.  This is only to make data generation the same on different operating systems: to fix an inconsistency with the mvrnorm function in the MASS package.
# group 1
# --------------------------------------
set.seed(200)
rs=.4
xy2 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
xy2 <- xy2[, order(apply(xy2, 2, var))[2:1]] #put the data columns in order from largest to smallest variance
# Note: Users do NOT have to re-order their data columns.  This is only to make data generation the same on different operating systems: to fix an inconsistency with the mvrnorm function in the MASS package.
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



Mx.expMean1Het <- matrix(c(0.03211284, -0.004889846),1,2)
Mx.expCov1Het <- matrix(c(1.0092856, 0.4813512, 0.4813512, 0.9935414),2,2)
Mx.expMean2Het <- matrix(c(0.03341992, -0.007112054),1,2)
Mx.expCov2Het <- matrix(c(1.012324, 0.3799160, 0.379916, 0.9956605),2,2)
Mx.LLHet <- 10944.873
# 1: Heterogeneity Model
# -------------------------------------

Mx.expMean1Hom <- matrix(c(0.03276872, -0.0059975),1,2)
Mx.expCov1Hom <- matrix(c(1.0108055, 0.4306339, 0.4306339, 0.9946009),2,2)
Mx.expMean2Hom <- matrix(c(0.03276872, -0.0059975),1,2)
Mx.expCov2Hom <- matrix(c(1.0108055, 0.4306339, 0.4306339, 0.9946009),2,2)
Mx.LLHom <- 10954.368
# 2: Homogeneity Model
# -------------------------------------
# Mx answers hard-coded
# -----------------------------------------------------------------------------


cov <- rbind(cbind(expCov1Het,expCov2Het),cbind(expCov1Hom,expCov2Hom))
mean <- rbind(cbind(expMean1Het, expMean2Het),cbind(expMean1Hom,expMean2Hom))
like <- rbind(cbind(LLHet),cbind(LLHom))
cov; mean; like
# OpenMx summary
# -----------------------------------------------------------------------------

Mx.cov <- rbind(cbind(Mx.expCov1Het,Mx.expCov2Het),cbind(Mx.expCov1Hom,Mx.expCov2Hom))
Mx.mean <- rbind(cbind(Mx.expMean1Het, Mx.expMean2Het),cbind(Mx.expMean1Hom,Mx.expMean2Hom))
Mx.like <- rbind(cbind(Mx.LLHet),cbind(Mx.LLHom))
Mx.cov; Mx.mean; Mx.like
# old Mx summary
# -----------------------------------------------------------------------------



omxCheckCloseEnough(LLHet,Mx.LLHet,.001)
omxCheckCloseEnough(expCov1Het,Mx.expCov1Het,.001)
omxCheckCloseEnough(expMean1Het,Mx.expMean1Het,.001)
omxCheckCloseEnough(expCov2Het,Mx.expCov2Het,.001)
omxCheckCloseEnough(expMean2Het,Mx.expMean2Het,.001)

omxCheckCloseEnough(LLHom,Mx.LLHom,.001)
omxCheckCloseEnough(expCov1Hom,Mx.expCov1Hom,.001)
omxCheckCloseEnough(expMean1Hom,Mx.expMean1Hom,.001)
omxCheckCloseEnough(expCov2Hom,Mx.expCov2Hom,.001)
omxCheckCloseEnough(expMean2Hom,Mx.expMean2Hom,.001)
# Compare OpenMx results to Mx results 
# (LL: likelihood; expCov: expected covariance, expMean: expected means)
# -----------------------------------------------------------------------------
