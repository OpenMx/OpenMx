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


#
# The virtual base class for all objective functions
#
setClass(Class = "MxBaseObjective", 
	representation = representation(
		name = "character",
		data = "MxCharOrNumber",
		dependencies = "character",
		result = "matrix", "VIRTUAL"))

setClassUnion("MxObjective", c("NULL", "MxBaseObjective"))

setGeneric("omxObjFunNamespace", 
	function(.Object, modelname, namespace) {
	return(standardGeneric("omxObjFunNamespace"))
})

setGeneric("omxObjFunConvert", 
	function(.Object, flatModel, model) {
	return(standardGeneric("omxObjFunConvert"))	
})

setGeneric("omxObjModelConvert",
	function(.Object, flatModel, model) {
	return(standardGeneric("omxObjModelConvert"))
})

setMethod("omxObjModelConvert", "MxBaseObjective",
	function(.Object, flatModel, model) {
		return(model)
})

setMethod("omxObjModelConvert", "NULL",
	function(.Object, flatModel, model) {
		return(model)
})

convertObjectives <- function(flatModel, model) {
	retval <- lapply(flatModel@objectives, function(x) {
		omxObjFunConvert(x, flatModel, model)
	})
	return(retval)
}

translateObjectives <- function(model, namespace) {
	flatModel <- omxFlattenModel(model, namespace)
	flatObjective <- flatModel@objectives[[omxIdentifier(model@name, 'objective')]]
	model <- omxObjModelConvert(flatObjective, flatModel, model)
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			submodel <- model@submodels[[i]]
			if (submodel@independent == FALSE) {
				model@submodels[[i]] <- translateObjectives(submodel, namespace)
			}
		}
	}
	return(model)
}
