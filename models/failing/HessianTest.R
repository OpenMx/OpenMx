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

# NOTE: The purpose of these tests is no longer to test the
# numDeriv derived Richardson extrapolation code. These tests
# are interesting because NPSOL and CSOLNP obtain different
# optimums.

library(OpenMx)

vNames <- c("V1", "V2", "V3", "V4")
dimNames <- list(2)
dimNames[[1]] <- vNames
dimNames[[2]] <- vNames
#dimNames
ObsLD  <- matrix(c( 1,  0,  0,  0,
                   .8,  1,  0,  0,
                   .7, .6,  1,  0,
                   .6, .5, .5,  1), nrow=4, byrow=TRUE)
ObsCov <- ObsLD %*% t(ObsLD)
#ObsCov
dimnames(ObsCov) <- dimNames
model <- mxModel("model",
                 mxMatrix(name="F", type="Full", free=TRUE, nrow=4, ncol=2),
                 mxMatrix(name="U", type="Diag", free=TRUE, nrow=4),
                 mxAlgebra(F %*% t(F) + U, name="PreCov", dimnames <- dimNames),
                 mxData(ObsCov, 'cov', numObs=500),
                 mxFitFunctionML(),mxExpectationNormal("PreCov"))
# start values
model$F@values <- matrix(c(.4, .5, .6, .2, .8, .7, .5, .2), nrow=4, ncol=2)
model$U@values <- diag(c(.2, .3, .4, .5))                            
# NOTE 10 observed statistics but 12 parameters. model is not identified                               
model <- mxOption(model, "Standard Errors", "Yes")
model <- mxRun(model)                
                                
# now examine the numerically differentiated hessian
fcn <- function(x) {
  F <- matrix(x[1:8], nrow=4, ncol=2)
  U <- diag(x[9:12])
  PreCov <- F %*% t(F) + U
  temp <- ObsCov * solve(PreCov)
  fval <- 499*(log(det(PreCov)) + sum(temp) - 4)
  return(fval)
}

data(numHess1)

omxCheckCloseEnough(data.matrix(numHess1), model@output$calculatedHessian, .01)

# use a different set of starting values
model2 <- mxModel("model2",
                 mxMatrix(name="F", type="Full", free=TRUE, nrow=4, ncol=2),
                 mxMatrix(name="U", type="Diag", free=TRUE, nrow=4),
                 mxAlgebra(F %*% t(F) + U, name="PreCov", dimnames <- dimNames),
                 mxData(ObsCov, 'cov', numObs=500),
                 mxFitFunctionML(),mxExpectationNormal("PreCov"))
model2$F@values <- matrix(c(.9, .8, .3, .1, .2, .4, .8, .9), nrow=4, ncol=2)
model2$U@values <- diag(c(.8, .6, .4, .9))
# again 10 obs stats & 12 parameters
model2 <- mxOption(model2, "Standard Errors", "Yes")
model2 <- mxRun(model2)

# numerical estimate of hessian
data(numHess2)

omxCheckCloseEnough(data.matrix(numHess2), model2@output$calculatedHessian, .01)
