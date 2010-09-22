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
# Program: LatentGrowthModel_PathRaw.R  
#  Author: Ryne Estabrook
#    Date: 09 17 2010 
#
# Growth Mixture Model
# Path style model input - Raw data input
#
# Revision History
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
data(myGrowthMixtureData)

#Create an MxModel object
# -----------------------------------------------------------------------
class1 <- mxModel("Class1", 
    mxMatrix(
        type="Full",
        nrow=7, 
        ncol=7,
        free=F,
        values=c(0,0,0,0,0,1,0,
                 0,0,0,0,0,1,1,
                 0,0,0,0,0,1,2,
                 0,0,0,0,0,1,3,
                 0,0,0,0,0,1,4,
                 0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0),
        byrow=TRUE,
        name="A"
    ),
    mxMatrix(
        type="Symm",
        nrow=7,
        ncol=7,
        free=c(T, F, F, F, F, F, F,
               F, T, F, F, F, F, F,
               F, F, T, F, F, F, F,
               F, F, F, T, F, F, F,
               F, F, F, F, T, F, F,
               F, F, F, F, F, T, T,
               F, F, F, F, F, T, T),
        values=c(0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  1,0.5,
                 0,0,0,0,0,0.5,  1),
        labels=c("residual", NA, NA, NA, NA, NA, NA,
                 NA, "residual", NA, NA, NA, NA, NA,
                 NA, NA, "residual", NA, NA, NA, NA,
                 NA, NA, NA, "residual", NA, NA, NA,
                 NA, NA, NA, NA, "residual", NA, NA,
                 NA, NA, NA, NA, NA, "vari1", "cov1",
                 NA, NA, NA, NA, NA, "cov1", "vars1"),
        byrow= TRUE,
        name="S"
    ),
    mxMatrix(
        type="Full",
        nrow=5,
        ncol=7,
        free=F,
        values=c(1,0,0,0,0,0,0,
                 0,1,0,0,0,0,0,
                 0,0,1,0,0,0,0,
                 0,0,0,1,0,0,0,
                 0,0,0,0,1,0,0),
        byrow=T,
        dimnames=list(NULL, c(names(myGrowthMixtureData), "intercept", "slope")),
        name="F"
    ),
    mxMatrix(
    	type="Full",
    	nrow=1, 
    	ncol=7,
        values=c(0,0,0,0,0,1,1),
        free=c(F,F,F,F,F,T,T),
        labels=c(NA,NA,NA,NA,NA,"meani1","means1"),
        name="M"
    ),
    mxRAMObjective("A","S","F","M", vector=TRUE)
) # close model

class2 <- mxModel(class1,
	mxMatrix(
        type="Symm",
        nrow=7,
        ncol=7,
        free=c(T, F, F, F, F, F, F,
               F, T, F, F, F, F, F,
               F, F, T, F, F, F, F,
               F, F, F, T, F, F, F,
               F, F, F, F, T, F, F,
               F, F, F, F, F, T, T,
               F, F, F, F, F, T, T),
        values=c(0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  0,  0,
                 0,0,0,0,0,  1,0.5,
                 0,0,0,0,0,0.5,  1),
        labels=c("residual", NA, NA, NA, NA, NA, NA,
                 NA, "residual", NA, NA, NA, NA, NA,
                 NA, NA, "residual", NA, NA, NA, NA,
                 NA, NA, NA, "residual", NA, NA, NA,
                 NA, NA, NA, NA, "residual", NA, NA,
                 NA, NA, NA, NA, NA, "vari2", "cov2",
                 NA, NA, NA, NA, NA, "cov2", "vars2"),
        byrow= TRUE,
        name="S"
    ),
    mxMatrix(
    	type="Full",
    	nrow=1, 
    	ncol=7,
        values=c(0,0,0,0,0,1,1),
        free=c(F,F,F,F,F,T,T),
        labels=c(NA,NA,NA,NA,NA,"meani2","means2"),
        name="M"
    ),
	name="Class2"
) # close model

     
# make a matrix of class probabilities
classP <- mxMatrix("Full", 2, 1, free=c(TRUE, FALSE), 
          values=.5, lbound=0.001, ubound=0.999,
          labels = c("pclass1", "pclass2[1,1]"), name="classProbs")

classA <- mxAlgebra(1-pclass1, name="pclass2")

algObj <- mxAlgebra(-2*sum(
          log(pclass1%x%Class1.objective + pclass2%x%Class2.objective)), 
          name="mixtureObj")

obj <- mxAlgebraObjective("mixtureObj")
      
gmm <- mxModel("Growth Mixture Model",
	mxData(
    	observed=myGrowthMixtureData,
        type="raw"
    ),
    class1, class2,
    classP, classA,
    algObj, obj
	)      

gmmFit <- mxRun(gmm)

summary(gmmFit)

gmmFit$classProbs

#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
#omxCheckCloseEnough(growthCurveFit@output$estimate[["meani"]], 9.930, 0.01)
