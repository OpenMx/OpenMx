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
	function(.Object, A, S, F, 
		data = as.numeric(NA), name = 'objective') {
		.Object@name <- name
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@data <- data
		return(.Object)
	}
)

setMethod("omxObjFunNamespace", signature("MxRAMObjective"), 
	function(.Object, modelname, namespace) {
		.Object@name <- omxIdentifier(modelname, .Object@name)
		.Object@A <- omxConvertIdentifier(.Object@A, modelname, namespace)
		.Object@S <- omxConvertIdentifier(.Object@S, modelname, namespace)
		.Object@F <- omxConvertIdentifier(.Object@F, modelname, namespace)
		.Object@data <- omxConvertIdentifier(.Object@data, modelname, namespace)
		return(.Object)
})

setMethod("omxObjFunConvert", signature("MxRAMObjective", "MxFlatModel"), 
	function(.Object, model, definitions) {
		name <- .Object@name
		aMatrix <- .Object@A
		sMatrix <- .Object@S
		fMatrix <- .Object@F
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The MxRAMObjective", omxQuotes(name),
				"does not have a dataset associated with it in model",
				omxQuotes(model@name))
			stop(msg, call. = FALSE)
		}
		.Object@A <- omxLocateIndex(model, aMatrix, name)
		.Object@S <- omxLocateIndex(model, sMatrix, name)
		.Object@F <- omxLocateIndex(model, fMatrix, name)
		.Object@data <- omxLocateIndex(model, data, name)
		return(.Object)
})

mxRAMObjective <- function(aMatrix = "A", sMatrix = "S", fMatrix = "F") {
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
	return(new("MxRAMObjective", aMatrix, sMatrix, fMatrix))
}
