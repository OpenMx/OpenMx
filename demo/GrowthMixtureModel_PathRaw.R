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
# Program: LatentGrowthModel_PathRaw.R  
# Author: Ryne Estabrook
# Date: 2010.09.17 
#
# ModelType: Growth Mixture
# DataType: Continuous
# Field: None
#
# Purpose: 
#      Growth Mixture Model
#      Path style model input - Raw data input
#
# RevisionHistory:
#      Ross Gore -- 2011.06.16 added Model, Data & Field metadata
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Libraries
# -----------------------------------------------------------------------------


data(myGrowthMixtureData)
# Prepare Data
# -----------------------------------------------------------------------------


class1 <- mxModel("Class1", 
    type="RAM",
    manifestVars=c("x1","x2","x3","x4","x5"),
    latentVars=c("intercept","slope"),
    # residual variances
    mxPath(
    	from=c("x1","x2","x3","x4","x5"), 
        arrows=2,
        free=TRUE, 
        values = c(1, 1, 1, 1, 1),
        labels=c("residual","residual","residual","residual","residual")
    ),
    # latent variances and covariance
    mxPath(
    	from=c("intercept","slope"), 
        arrows=2,
		connect="unique.pairs",
        free=TRUE, 
        values=c(1, .4, 1),
        labels=c("vari1", "cov1", "vars1")
    ),
    # intercept loadings
    mxPath(
    	from="intercept",
        to=c("x1","x2","x3","x4","x5"),
        arrows=1,
        free=FALSE,
        values=c(1, 1, 1, 1, 1)
    ),
    # slope loadings
    mxPath(
    	from="slope",
        to=c("x1","x2","x3","x4","x5"),
        arrows=1,
        free=FALSE,
        values=c(0, 1, 2, 3, 4)
    ),
    # manifest means
    mxPath(from="one",
        to=c("x1", "x2", "x3", "x4", "x5"),
        arrows=1,
        free=FALSE,
        values=c(0, 0, 0, 0, 0)
    ),
    # latent means
    mxPath(from="one",
        to=c("intercept", "slope"),
        arrows=1,
        free=TRUE,
        values=c(0, -1),
        labels=c("meani1", "means1")
    ),
    # enable the likelihood vector
    mxRAMObjective(A = "A",
        S = "S",
        F = "F",
        M = "M",
        vector = TRUE)
) # close model

class2 <- mxModel(class1,
	# latent variances and covariance
    mxPath(
    	from=c("intercept","slope"), 
        arrows=2,
		connect="unique.pairs",
        free=TRUE, 
        values=c(1, .5, 1),
        labels=c("vari2", "cov2", "vars2")
    ),
    # latent means
    mxPath(from="one",
        to=c("intercept", "slope"),
        arrows=1,
        free=TRUE,
        values=c(5, 1),
        labels=c("meani2", "means2")
    ),
	name="Class2"
) # close model
# Create an MxModel object
# -----------------------------------------------------------------------------

# request that individual likelihoods are used
# required for correct parameterization of class probabilities
# -----------------------------------------------------------------------------     

classP <- mxMatrix("Full", 2, 1, free=c(TRUE, FALSE), 
          values=1, lbound=0.001, 
          labels = c("p1", "p2"), name="Props")
# make a matrix of class probabilities
# -----------------------------------------------------------------------------

classS <- mxAlgebra(Props%x%(1/sum(Props)), name="classProbs")

algObj <- mxAlgebra(-2*sum(
          log(classProbs[1,1]%x%Class1.objective + classProbs[2,1]%x%Class2.objective)), 
          name="mixtureObj")
          
obj <- mxAlgebraObjective("mixtureObj")
      
gmm <- mxModel("Growth Mixture Model",
	mxData(
    	observed=myGrowthMixtureData,
        type="raw"
    ),
    class1, class2,
    classP, classS,
    algObj, obj
	)      

gmmFit <- mxRun(gmm, suppressWarnings=TRUE)

summary(gmmFit)

gmmFit$classProbs


omxCheckCloseEnough(gmmFit@output$Minus2LogLikelihood, 8739.05, 0.01)
omxCheckCloseEnough(max(mxEval(classProbs, gmmFit)), 0.6009, 0.01)
omxCheckCloseEnough(min(mxEval(classProbs, gmmFit)), 0.3991, 0.01)
# Check to see if results match within the specified bounds
# -----------------------------------------------------------------------------
