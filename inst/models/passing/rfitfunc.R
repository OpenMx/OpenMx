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


library(OpenMx)

cubic <- function(marg, state){
    x <- marg$matrices$param$values[1]
    got <- (x-5)*(x-1)*(x+3)
    return(got)
}

model <- mxModel(name="root",
	      mxMatrix(type="Full", ncol=1, nrow=1, name="param", free=TRUE, values=0),
	      mxFitFunctionR(cubic))
model <- mxRun(model, silent=TRUE, suppressWarnings=TRUE)
omxCheckCloseEnough(model$matrices$param$values, 3.309401, 10^-3)

###

infer <- function(marg,state) return(Inf)

model <- mxModel(name="inf",
		 mxMatrix(type="Full", ncol=1, nrow=1, name="param", free=TRUE, values=0),
		 mxFitFunctionR(infer))
ign <- omxCheckWarning(try(mxRun(model), silent=TRUE),
		"In model 'inf' Optimizer returned a non-zero status code 10. Starting values are not feasible. Consider mxTryHard()")

###

counter <<- 1
count <- function(marg,state) {
	counter <<- state[[1]]
	state[[1]] <- state[[1]] + 1
	return(list(1, state))
}

model <- mxModel(name="count",
		 mxMatrix(type="Full", ncol=1, nrow=1, name="param", free=TRUE, values=0),
		 mxFitFunctionR(count, 1))
model <- omxCheckWarning(mxRun(model, silent=TRUE),
                         "In model 'count' Optimizer returned a non-zero status code 5. The Hessian at the solution does not appear to be convex. See ?mxCheckIdentification for possible diagnosis (Mx status RED).")

omxCheckTrue(counter > 1)

###

omxCheckError(mxCheckIdentification(model), "Identification check is not possible for models with 'MxFitFunctionAlgebra', 'MxFitFunctionRow', and 'MxFitFunctionR' fit functions.
 If you have a multigroup model, use mxFitFunctionMultigroup.")


###

toomany <- function(marg,state) {
	return(list(1, state, 5))
}

model <- mxModel(name="toomany",
		 mxMatrix(type="Full", ncol=1, nrow=1, name="param", free=TRUE, values=0),
		 mxFitFunctionR(toomany))
model <- omxCheckWarning(
  omxCheckError(mxRun(model), "The job for model 'toomany' exited abnormally with the error message: FitFunction returned more than 2 arguments"),
  "In model 'toomany' Optimizer returned a non-zero status code 10. Starting values are not feasible. Consider mxTryHard()")
