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


setClass(Class = "MxRAMObjective",
	representation = representation(
		A = "MxCharOrNumber",
		S = "MxCharOrNumber",
		F = "MxCharOrNumber"),
	contains = "MxBaseObjective")

setMethod("initialize", "MxRAMObjective",
	function(.Object, name, A, S, F, data = NA_real_) {
		.Object@name <- name
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@data <- data
		return(.Object)
	}
)

setMethod("omxObjFunConvert", signature("MxRAMObjective", "MxFlatModel"), 
	function(.Object, model, definitions) {
		name <- .Object@name
		aMatrix <- .Object@A
		sMatrix <- .Object@S
		fMatrix <- .Object@F
		data <- .Object@data
		A <- omxLocateIndex(model, aMatrix, name)
		S <- omxLocateIndex(model, sMatrix, name)
		F <- omxLocateIndex(model, fMatrix, name)
		dIndex <- omxLocateIndex(model, data, name)
		return(new("MxRAMObjective", name, A, S, F, dIndex))
})

mxRAMObjective <- function(aMatrix = "A", sMatrix = "S", fMatrix = "F", name = NA) {
	if(is.na(name)) name <- omxUntitledName()
	if (typeof(name) != "character") {
		msg <- paste("Name argument is not a string",
			"(the name of the objective function)")
		stop(msg)
	}
	if (typeof(aMatrix) != "character") {
		msg <- paste("aMatrix argument is not a string",
			"(the name of the 'A' matrix)")
		stop(msg)
	}	
	if (typeof(sMatrix) != "character") {
		msg <- paste("sMatrix argument is not a string",
			"(the name of the 'S' matrix)")
		stop(msg)
	}
	if (typeof(fMatrix) != "character") {
		msg <- paste("fMatrix argument is not a string",
			"(the name of the 'F' matrix)")
		stop(msg)
	}
	return(new("MxRAMObjective", name, aMatrix, sMatrix, fMatrix))
}
