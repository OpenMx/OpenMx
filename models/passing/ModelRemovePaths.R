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


# -----------------------------------------------------------------------
# Program: TwoFactorModel_PathCov.R  
#  Author: Ryne Estabrook
#    Date: 08 01 2009 
#
# Two Factor model to estimate factor loadings, residual variances and means
# Path style model input - Covariance matrix data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
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

#Create an MxModel object
# -----------------------------------------------------------------------
twoFactorModel <- mxModel("Two Factor Model Path", type="RAM",
    mxData(
    	observed=twoFactorCov, 
    	type="cov", 
    	numObs=500, 
    	means=twoFactorMeans
    ),
    manifestVars=c("x1", "x2", "x3", "y1", "y2", "y3"),
    latentVars=c("F1","F2"),
    # residual variances
    mxPath(
    	from=c("x1", "x2", "x3", "y1", "y2", "y3"),
        arrows=2, 
        free=TRUE, 
        values=1,
        labels=c("e1","e2","e3","e4","e5","e6")
    ),
    # latent variances and covaraince
    mxPath(
    	from=c("F1","F2"),
        arrows=2,
#        all=TRUE, 
        connect="unique.pairs",
        free=TRUE,
        values=c(1, .5, 1),
        labels=c("varF1", "cov", "varF2")
    ), 
    # factor loadings for x variables
    mxPath(
    	from="F1",
        to=c("x1","x2","x3"),
        arrows=1,
        free=c(FALSE,TRUE,TRUE),
        values=c(1,1,1),
        labels=c("l1","l2","l3")
    ),
    # factor loadings for y variables
    mxPath(
    	from="F2",
        to=c("y1","y2","y3"),
        arrows=1,
        free=c(FALSE,TRUE,TRUE),
        values=c(1,1,1),
        labels=c("l4","l5","l6")
    ),
    # means
    mxPath(
    	from="one",
        to=c("x1","x2","x3","y1","y2","y3","F1","F2"),
        arrows=1,
        free=c(T ,T, T, T, T, T, F, F),
        values=c(1,1,1,1,1,1,0,0),
        labels=c("meanx1","meanx2","meanx3","meany1","meany2","meany3",NA,NA)
    )
) # close model

elimination <- twoFactorModel

elimination <- mxModel(elimination, remove = TRUE,
 	mxPath(
        from=c("x1", "x2", "x3", "y1", "y2", "y3"),
		arrows=2
    ),
    mxPath(
        from=c("F1","F2"),
        arrows=2,
#        all=TRUE
        connect="unique.pairs"
    ),
	mxPath(
        from="F1",
        to=c("x1","x2","x3"),
        arrows=1
    ),
    mxPath(
        from="F2",
        to=c("y1","y2","y3"),
        arrows=1
    ),
    # means
    mxPath(
        from="one",
        to=c("x1","x2","x3","y1","y2","y3","F1","F2")
    )
)

# Don't ever do this to your scripts.
# We need to strip the dimnames so that the
# checking below can ignore them.
dimnames(elimination$A$values) <- NULL
dimnames(elimination$A$free) <- NULL
dimnames(elimination$A$labels) <- NULL

dimnames(elimination$S$values) <- NULL
dimnames(elimination$S$free) <- NULL
dimnames(elimination$S$labels) <- NULL

dimnames(elimination$M$values) <- NULL
dimnames(elimination$M$free) <- NULL
dimnames(elimination$M$labels) <- NULL

omxCheckIdentical(elimination$A$values, matrix(0, 8, 8))
omxCheckIdentical(elimination$A$free, matrix(FALSE, 8, 8))
omxCheckIdentical(elimination$A$labels, matrix(as.character(NA), 8, 8))

omxCheckIdentical(elimination$S$values, matrix(0, 8, 8))
omxCheckIdentical(elimination$S$free, matrix(FALSE, 8, 8))
omxCheckIdentical(elimination$S$labels, matrix(as.character(NA), 8, 8))

omxCheckIdentical(elimination$M$values, matrix(0, 1, 8))
omxCheckIdentical(elimination$M$free, matrix(FALSE, 1, 8))
omxCheckIdentical(elimination$M$labels, matrix(as.character(NA), 1, 8))
