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


# -----------------------------------------------------------------------------
# Program: OneFactorModel_PathRaw.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Factor
# DataType: Continuous
# Field: None
#
# Purpose: 
#      One Factor model to estimate factor loadings, residual variances and means
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Hermine Maes -- 2009.10.08 updated & reformatted
#      Ross Gore -- 2011.06.06 added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

data(myFADataRaw)
# Prepare Data
# -----------------------------------------------------------------------------

myFADataRaw <- myFADataRaw[,c("x1","x2","x3","x4","x5","x6")]

oneFactorModel <- mxModel("parent",
	mxAlgebra(-2 * sum(log(likelihoods.objective)), name = "algObjective"),
	mxAlgebraObjective("algObjective"),
	mxModel("likelihoods", 
		type="RAM",
		mxData(
			observed=myFADataRaw, 
			type="raw"
		),
		manifestVars=c("x1","x2","x3","x4","x5","x6"),
		latentVars="F1",
		mxPath(from=c("x1","x2","x3","x4","x5","x6"),
			arrows=2,
			free=TRUE,
			values=c(1,1,1,1,1,1),
			labels=c("e1","e2","e3","e4","e5","e6")
		), 
		# residual variances
		# -------------------------------------
		mxPath(from="F1",
			arrows=2,
			free=TRUE,
			values=1,
			labels ="varF1"
		), 
		# latent variance
		# -------------------------------------
		mxPath(from="F1",
			to=c("x1","x2","x3","x4","x5","x6"),
			arrows=1,
			free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
			values=c(1,1,1,1,1,1),
			labels =c("l1","l2","l3","l4","l5","l6")
		), 
		# factor loadings
		# -------------------------------------
		mxPath(from="one",
			to=c("x1","x2","x3","x4","x5","x6","F1"),
			arrows=1,
			free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE),
			values=c(1,1,1,1,1,1,0),
			labels =c("meanx1","meanx2","meanx3","meanx4","meanx5","meanx6",NA)
		),
		# means
		# -------------------------------------
		mxRAMObjective("A", "S", "F", "M", vector=TRUE)
	) # close model 'likelihoods'
) #close model 'parent'
# Create an MxModel object
# -----------------------------------------------------------------------------
   
oneFactorFit <- mxRun(oneFactorModel)      

summary(oneFactorFit)
oneFactorFit@output$estimate


omxCheckCloseEnough(oneFactorFit@output$estimate[["l2"]], 0.999, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l3"]], 0.959, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l4"]], 1.028, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l5"]], 1.008, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l6"]], 1.021, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["varF1"]], 0.645, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e1"]], 0.350, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e2"]], 0.379, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e3"]], 0.389, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e4"]], 0.320, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e5"]], 0.370, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e6"]], 0.346, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx1"]], 2.988, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx2"]], 3.011, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx3"]], 2.986, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx4"]], 3.053, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx5"]], 3.016, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx6"]], 3.010, 0.01)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------