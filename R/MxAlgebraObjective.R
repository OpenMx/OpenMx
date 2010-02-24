#
#   Copyright 2007-2009 The OpenMx Project
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


setClass(Class = "MxAlgebraObjective",
	representation = representation(
		algebra = "MxCharOrNumber"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxAlgebraObjective",
	function(.Object, algebra, 
		data = as.integer(NA), name = 'objective') {
		.Object@name <- name
		.Object@algebra <- algebra
		.Object@data <- data
		return(.Object)
	}
)

setMethod("genericObjDependencies", signature("MxAlgebraObjective"),
	function(.Object, dependencies) {
	dependencies <- omxAddDependency(.Object@algebra, .Object@name, dependencies)
	return(dependencies)
})

setMethod("genericObjFunConvert", signature("MxAlgebraObjective"), 
	function(.Object, flatModel, model, defVars) {
		name <- .Object@name
		algebra <- .Object@algebra
		if (is.na(algebra)) {
			modelname <- omxReverseIdentifier(model, .Object@name)[[1]]
			msg <- paste("The algebra name cannot be NA",
			"for the algebra objective of model", omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		algebraIndex <- omxLocateIndex(flatModel, algebra, name)
		.Object@algebra <- algebraIndex
		return(.Object)
})

setMethod("genericObjFunNamespace", signature("MxAlgebraObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		.Object@algebra <- omxConvertIdentifier(.Object@algebra, modelname, namespace)
		return(.Object)
})

setMethod("genericObjRename", signature("MxAlgebraObjective"),
	function(.Object, oldname, newname) {
		.Object@algebra <- renameReference(.Object@algebra, oldname, newname)
		return(.Object)
})

mxAlgebraObjective <- function(algebra) {
	if (missing(algebra) || typeof(algebra) != "character") {
		stop("Algebra argument is not a string (the name of the algebra)")
	}
	return(new("MxAlgebraObjective", algebra))
}

displayAlgebraObjective <- function(objective) {
	cat("MxAlgebraObjective", omxQuotes(objective@name), '\n')
	cat("@algebra\n", omxQuotes(objective@algebra), '\n')
	if (length(objective@result) == 0) {
		cat("@result: (not yet computed) ")
	} else {
		cat("@result:\n")
	}
	print(objective@result)
	invisible(objective)
}

setMethod("print", "MxAlgebraObjective", function(x,...) { 
	displayAlgebraObjective(x) 
})

setMethod("show", "MxAlgebraObjective", function(object) { 
	displayAlgebraObjective(object) 
})
