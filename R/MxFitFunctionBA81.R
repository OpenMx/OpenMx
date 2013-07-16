#
#   Copyright 2012 The OpenMx Project
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


setClass(Class = "MxFitFunctionBA81",
	representation = representation(
          data = "MxCharOrNumber"),
	contains = "MxBaseFitFunction")

setMethod("initialize", "MxFitFunctionBA81",
	function(.Object) {
                .Object@name <- 'fitfunction'
		.Object@data <- as.integer(NA)
		return(.Object)
	}
)

setMethod("genericFitFunConvert", signature("MxFitFunctionBA81"),
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		name <- .Object@name
		modelname <- imxReverseIdentifier(model, name)[[1]]
		expectName <- paste(modelname, "expectation", sep=".")
		if (expectName %in% names(flatModel@expectations)) {
			expectIndex <- imxLocateIndex(flatModel, expectName, name)
		} else {
			expectIndex <- as.integer(NA)
		}
		.Object@expectation <- expectIndex

                for (s in c("data")) {
			if (is.null(slot(.Object, s))) next
			slot(.Object, s) <-
			  imxLocateIndex(flatModel, slot(.Object, s), name)
                }
		return(.Object)
})

setMethod("qualifyNames", signature("MxFitFunctionBA81"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		return(.Object)
})

setMethod("genericFitRename", signature("MxFitFunctionBA81"),
	function(.Object, oldname, newname) {
          # not sure what goes here yet
		return(.Object)
})

mxFitFunctionBA81 <- function() {
	return(new("MxFitFunctionBA81"))
}
