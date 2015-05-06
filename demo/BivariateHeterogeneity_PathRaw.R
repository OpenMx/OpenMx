#
#   Copyright 2007-2015 The OpenMx Project
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
# Program: BivariateHeterogeneity_PathRaw.R  
# Author: Hermine Maes
# Date: 2009.08.01 
#
# ModelType: Heterogeneity
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Bivariate Heterogeneity model to test differences in means and variances across multiple groups
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field metadata
#      Hermine Maes -- 2014.11.02 piecewise specification
# -----------------------------------------------------------------------

require(OpenMx)
require(MASS)
# Load Libraries
# -----------------------------------------------------------------------------


if (0) {
	set.seed(200)
	rs=.5
	xy1 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
	set.seed(200)
	# group 1
	# -------------------------------------

	rs=.4
	xy2 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
	# group 2
	# -------------------------------------

	xy1 <- round(xy1, 12)
	write.csv(xy1, "data/BivariateHeterogeneity1.csv", row.names=FALSE)
	xy2 <- round(xy2, 12)
	write.csv(xy2, "data/BivariateHeterogeneity2.csv", row.names=FALSE)
	stop("data generated")
} else {
	xy1 <- as.matrix(read.csv("data/BivariateHeterogeneity1.csv"))
	xy2 <- as.matrix(read.csv("data/BivariateHeterogeneity2.csv"))
}

selVars <- c("X","Y")
summary(xy1)
cov(xy1)
dimnames(xy1) <- list(NULL, selVars)
summary(xy2)
cov(xy2)
dimnames(xy2) <- list(NULL, selVars)
# Print Descriptive Statistics
# -------------------------------------
# Simulate Data
# -----------------------------------------------------------------------------

dataRaw1     <- mxData( observed=xy1, type="raw")
variances1   <- mxPath( from=selVars, arrows=2, 
          	            free=T, values=1, lbound=.01, labels=c("vX1","vY1") )
covariance1  <- mxPath( from="X", to="Y", arrows=2, 
            	        free=T, values=.2, lbound=.01, labels="cXY1")
means1       <- mxPath( from="one", to=selVars, arrows=1, 
        		        free=T, values=c(0.1,-0.1), ubound=c(NA,0), lbound=c(0,NA), 
        		        labels=c("mX1","mY1") )
model1       <- mxModel("group1", type="RAM", manifestVars= selVars,
                         dataRaw1, variances1, covariance1, means1)

dataRaw2     <- mxData( observed=xy2, type="raw")
variances2   <- mxPath( from=selVars, arrows=2, 
          	            free=T, values=1, lbound=.01, labels=c("vX2","vY2") )
covariance2  <- mxPath( from="X", to="Y", arrows=2, 
            	        free=T, values=.2, lbound=.01, labels="cXY2")
means2       <- mxPath( from="one", to=selVars, arrows=1, 
        		        free=T, values=c(0.1,-0.1), ubound=c(NA,0), lbound=c(0,NA), 
        		        labels=c("mX2","mY2") )
model2       <- mxModel("group2", type="RAM", manifestVars= selVars,
                         dataRaw2, variances2, covariance2, means2)

#h12          <- mxAlgebra( group1.fitfunction + group2.fitfunction, name="h12" )
#funML        <- mxFitFunctionAlgebra("h12")
fun           <- mxFitFunctionMultigroup(c("group1.fitfunction", "group2.fitfunction"))
bivHetModel   <- mxModel("bivariate Heterogeneity Path Specification",
                        model1, model2, fun )

bivHetFit <- mxRun(bivHetModel)
EM1Het <- bivHetFit$group1.fitfunction$info$expMean
EM2Het <- bivHetFit$group2.fitfunction$info$expMean
EC1Het <- bivHetFit$group1.fitfunction$info$expCov
EC2Het <- bivHetFit$group2.fitfunction$info$expCov
LLHet <- bivHetFit$output$fit

EM1Het; EM2Het; EC1Het; EC2Het; LLHet

# Fit Heterogeneity Model
# -----------------------------------------------------------------------------

bivHomModel <- bivHetModel
bivHomModel[['group2.S']]$labels <- bivHomModel[['group1.S']]$labels
bivHomModel[['group2.M']]$labels <- bivHomModel[['group1.M']]$labels

bivHomFit <- mxRun(bivHomModel)
EM1Hom <- bivHomFit$group1.fitfunction$info$expMean
EM2Hom <- bivHomFit$group2.fitfunction$info$expMean
EC1Hom <- bivHomFit$group1.fitfunction$info$expCov
EC2Hom <- bivHomFit$group2.fitfunction$info$expCov
LLHom <- bivHomFit$output$fit

EM1Hom; EM2Hom; EC1Hom; EC2Hom; LLHom

Chi= LLHom-LLHet
LRT= rbind(LLHet,LLHom,Chi)
LRT
# Fit Homnogeneity Model
# -----------------------------------------------------------------------------



Mx.EM1Het <- matrix(c(0.03211284, -0.004889846),1,2)
Mx.EC1Het <- matrix(c(1.0092856, 0.4813512, 0.4813512, 0.9935414),2,2)
Mx.EM2Het <- matrix(c(0.03341992, -0.007112054),1,2)
Mx.EC2Het <- matrix(c(1.012324, 0.3799160, 0.379916, 0.9956605),2,2)
Mx.LLHet <- 10944.873
# 1: Heterogeneity Model
# -------------------------------------
Mx.EM1Hom <- matrix(c(0.03276872, -0.0059975),1,2)
Mx.EC1Hom <- matrix(c(1.0108055, 0.4306339, 0.4306339, 0.9946009),2,2)
Mx.EM2Hom <- matrix(c(0.03276872, -0.0059975),1,2)
Mx.EC2Hom <- matrix(c(1.0108055, 0.4306339, 0.4306339, 0.9946009),2,2)
Mx.LLHom <- 10954.368
# 2: Homogeneity Model
# -------------------------------------
# Mx answers hard-coded
# -----------------------------------------------------------------------------

cov <- rbind(cbind(EC1Het,EC2Het),cbind(EC1Hom,EC2Hom))
mean <- rbind(cbind(EM1Het, EM2Het),cbind(EM1Hom,EM2Hom))
like <- rbind(cbind(LLHet),cbind(LLHom))
cov; mean; like
# OpenMx summary
# -----------------------------------------------------------------------------

Mx.cov <- rbind(cbind(Mx.EC1Het,Mx.EC2Het),cbind(Mx.EC1Hom,Mx.EC2Hom))
Mx.mean <- rbind(cbind(Mx.EM1Het, Mx.EM2Het),cbind(Mx.EM1Hom,Mx.EM2Hom))
Mx.like <- rbind(cbind(Mx.LLHet),cbind(Mx.LLHom))
Mx.cov; Mx.mean; Mx.like
# old Mx summary
# -----------------------------------------------------------------------------


omxCheckCloseEnough(LLHet,Mx.LLHet,.01)
omxCheckCloseEnough(EC1Het,Mx.EC1Het,.001)
omxCheckCloseEnough(EM1Het,Mx.EM1Het,.001)
omxCheckCloseEnough(EC2Het,Mx.EC2Het,.001)
omxCheckCloseEnough(EM2Het,Mx.EM2Het,.001)

omxCheckCloseEnough(LLHom,Mx.LLHom,.001)
omxCheckCloseEnough(EC1Hom,Mx.EC1Hom,.001)
omxCheckCloseEnough(EM1Hom,Mx.EM1Hom,.001)
omxCheckCloseEnough(EC2Hom,Mx.EC2Hom,.001)
omxCheckCloseEnough(EM2Hom,Mx.EM2Hom,.001)
# Compare OpenMx results to Mx results 
# (LL: likelihood; EC: expected covariance, EM: expected means)
# -----------------------------------------------------------------------------
