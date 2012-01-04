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


# -----------------------------------------------------------------------
# Program: TwoFactorModel_PathCov.R  
# Author: Ryne Estabrook
# Date: 2009.08.01 
#
# ModelType: Factor
# DataType: Continuous
# Field: None
#
# Purpose: 
# 	Two Factor model to estimate factor loadings, residual variances and means
# 	Path style model input - Covariance matrix data input
#
# RevisionHistory:
# 	Hermine Maes -- 2009.10.08 updated & reformatted
# 	Ross Gore -- 2011.06.06	added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------

myFADataCov <- matrix(
      c(0.997, 0.642, 0.611, 0.672, 0.637, 0.677, 0.342, 0.299, 0.337,
        0.642, 1.025, 0.608, 0.668, 0.643, 0.676, 0.273, 0.282, 0.287,
        0.611, 0.608, 0.984, 0.633, 0.657, 0.626, 0.286, 0.287, 0.264,
        0.672, 0.668, 0.633, 1.003, 0.676, 0.665, 0.330, 0.290, 0.274,
        0.637, 0.643, 0.657, 0.676, 1.028, 0.654, 0.328, 0.317, 0.331,
        0.677, 0.676, 0.626, 0.665, 0.654, 1.020, 0.323, 0.341, 0.349,
        0.342, 0.273, 0.286, 0.330, 0.328, 0.323, 0.993, 0.472, 0.467,
        0.299, 0.282, 0.287, 0.290, 0.317, 0.341, 0.472, 0.978, 0.507,
        0.337, 0.287, 0.264, 0.274, 0.331, 0.349, 0.467, 0.507, 1.059),
      nrow=9,
      dimnames=list(
          c("x1", "x2", "x3", "x4", "x5", "x6", "y1", "y2", "y3"),
          c("x1", "x2", "x3", "x4", "x5", "x6", "y1", "y2", "y3")),
)

twoFactorCov <- myFADataCov[c("x1","x2","x3","y1","y2","y3"),c("x1","x2","x3","y1","y2","y3")]
  
myFADataMeans <- c(2.988, 3.011, 2.986, 3.053, 3.016, 3.010, 2.955, 2.956, 2.967)
names(myFADataMeans) <- c("x1", "x2", "x3", "x4", "x5", "x6", "y1", "y2", "y3")
  
twoFactorMeans <- myFADataMeans[c(1:3,7:9)]
# Prepare Data
# -----------------------------------------------------------------------------

twoFactorModel <- mxModel("Two Factor Model Path", type="RAM",
    mxData(
    	observed=twoFactorCov, 
    	type="cov", 
    	numObs=500, 
    	means=twoFactorMeans
    ),
    manifestVars=c("x1", "x2", "x3", "y1", "y2", "y3"),
    latentVars=c("F1","F2"),
    mxPath(
    	from=c("x1", "x2", "x3", "y1", "y2", "y3"),
        arrows=2, 
        free=TRUE, 
        values=1,
        labels=c("e1","e2","e3","e4","e5","e6")
    ),
	# residual variances
	# -------------------------------------
    mxPath(
    	from=c("F1","F2"),
        arrows=2,
		connect="unique.pairs",
        free=TRUE,
        values=c(1, .5, 1),
        labels=c("varF1", "cov", "varF2")
    ), 
	# latent variances and covaraince
	# -------------------------------------
    mxPath(
    	from="F1",
        to=c("x1","x2","x3"),
        arrows=1,
        free=c(FALSE,TRUE,TRUE),
        values=c(1,1,1),
        labels=c("l1","l2","l3")
    ),
	# factor loadings for x variables
	# -------------------------------------
    mxPath(
    	from="F2",
        to=c("y1","y2","y3"),
        arrows=1,
        free=c(FALSE,TRUE,TRUE),
        values=c(1,1,1),
        labels=c("l4","l5","l6")
    ),
	# factor loadings for y variables
	# -------------------------------------
    mxPath(
    	from="one",
        to=c("x1","x2","x3","y1","y2","y3","F1","F2"),
        arrows=1,
        free=c(T ,T, T, T, T, T, F, F),
        values=c(1,1,1,1,1,1,0,0),
        labels=c("meanx1","meanx2","meanx3","meany1","meany2","meany3",NA,NA)
    )
	# means
	# -------------------------------------
) # close model
# Create an MxModel object
# -----------------------------------------------------------------------------
      
twoFactorFit <- mxRun(twoFactorModel)

summary(twoFactorFit)
twoFactorFit@output$estimate

omxCheckCloseEnough(twoFactorFit@output$estimate[["l2"]], 0.9720, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l3"]], 0.9310, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l5"]], 1.0498, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l6"]], 1.0533, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["varF1"]], 0.6622, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["varF2"]], 0.4510, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["cov"]], 0.2958, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e1"]], 0.3348, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e2"]], 0.3994, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e3"]], 0.4101, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e4"]], 0.5420, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e5"]], 0.4809, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e6"]], 0.5586, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx1"]], 2.988, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx2"]], 3.011, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx3"]], 2.986, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meany1"]], 2.955, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meany2"]], 2.956, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meany3"]], 2.967, 0.01)
# Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------------