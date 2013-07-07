#
#   Copyright 2013 The OpenMx Project
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

setClass(Class = "MxFitFunctionMultigroup",
	 representation = representation(
	   groups = "MxOptionalCharOrNumber"),
	 contains = "MxBaseFitFunction")

setMethod("initialize", "MxFitFunctionMultigroup",
	function(.Object, groups, name = 'fitfunction') {
		.Object@name <- name
		.Object@groups <- groups
		return(.Object)
	}
)

setMethod("genericFitDependencies", signature("MxFitFunctionMultigroup"),
	function(.Object, flatModel, dependencies) {
	dependencies <- callNextMethod()
	dependencies <- imxAddDependency(.Object@groups, .Object@name, dependencies)
	return(dependencies)
})

setMethod("qualifyNames", signature("MxFitFunctionMultigroup"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		return(.Object)
})

# "model.algebra" or "model" for "model.fitfunction"
setMethod("genericFitFunConvert", "MxFitFunctionMultigroup", 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
		name <- .Object@name
		if (length(.Object@groups)) .Object@groups <- vapply(.Object@groups, function(group) {
			path <- unlist(strsplit(group, imxSeparatorChar, fixed = TRUE))
			if (length(path) == 1) {
				ff <- paste(path, "fitfunction", sep=".")
				length(model@algebras) + imxLocateIndex(flatModel, ff, name)
			} else if (length(path) == 2) {
				# restrict to algebra or fitfunction TODO
				imxLocateIndex(flatModel, group, name)
			}
		}, 1L)
		return(.Object)
})


mxFitFunctionMultigroup <- function(groups) {
	return(new("MxFitFunctionMultigroup", groups))
}
