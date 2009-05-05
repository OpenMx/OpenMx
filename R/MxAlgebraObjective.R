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
	function(.Object, name, algebra) {
		.Object@name <- name
		.Object@algebra <- algebra
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxAlgebraObjective", "MxFlatModel"), function(.Object, model, definitions) {
		name <- .Object@name
		algebra <- .Object@algebra
		algebraIndex <- omxLocateIndex(model, algebra, name)
		if (is.na(algebraIndex)) {
			stop(paste("Could not find a matrix/algebra with name", 
				algebra, "in the model."))
		}
		return(new("MxAlgebraObjective", name, algebraIndex))
})


mxAlgebraObjective <- function(algebra, name = NA) {
	if(is.na(name)) name <- omxUntitledName()
	omxVerifyName(name)
	if (typeof(name) != "character") {
		stop("Name argument is not a string (the name of the objective function)")
	}
	if (missing(algebra) || typeof(algebra) != "character") {
		stop("Algebra argument is not a string (the name of the algebra)")
	}
	return(new("MxAlgebraObjective", name, algebra))
}
