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


#
# The virtual base class for all objective functions
#
setClass(Class = "MxBaseObjective", 
	representation = representation(
		name = "character",
		data = "MxCharOrNumber",
		result = "matrix", "VIRTUAL"))

setClassUnion("MxObjective", c("NULL", "MxBaseObjective"))

setGeneric("genericObjFunNamespace", 
	function(.Object, modelname, namespace) {
	return(standardGeneric("genericObjFunNamespace"))
})

setGeneric("genericObjFunConvert", 
	function(.Object, flatModel, model) {
	return(standardGeneric("genericObjFunConvert"))	
})

setGeneric("genericObjModelConvert",
	function(.Object, job, model, flatJob) {
	return(standardGeneric("genericObjModelConvert"))
})

setGeneric("genericObjDependencies",
	function(.Object, dependencies) {
	return(standardGeneric("genericObjDependencies"))
})

setGeneric("genericObjInitialMatrix",
	function(.Object, flatModel) {
	return(standardGeneric("genericObjInitialMatrix"))
})

setGeneric("genericObjNewEntities",
	function(.Object) {
	return(standardGeneric("genericObjNewEntities"))
})

setGeneric("genericObjRename",
	function(.Object, oldname, newname) {
	return(standardGeneric("genericObjRename"))
})

setMethod("genericObjModelConvert", "MxBaseObjective",
	function(.Object, job, model, flatJob) {
		return(job)
})

setMethod("genericObjModelConvert", "NULL",
	function(.Object, job, model, flatJob) {
		return(job)
})

setMethod("genericObjDependencies", "MxBaseObjective",
	function(.Object, dependencies) {
		return(dependencies)
})

setMethod("genericObjDependencies", "NULL",
	function(.Object, dependencies) {
		return(dependencies)
})

setMethod("genericObjNewEntities", "MxBaseObjective",
	function(.Object) {
		return(NULL)
})

setMethod("genericObjInitialMatrix", "MxBaseObjective",
	function(.Object, flatModel) {
		return(matrix(as.double(NA), 1, 1))
})

setMethod("genericObjInitialMatrix", "NULL",
	function(.Object, flatModel) {
		return(NULL)
})

setMethod("genericObjRename", "MxBaseObjective",
	function(.Object, oldname, newname) {
		return(.Object)
})

setMethod("genericObjRename", "NULL",
	function(.Object, oldname, newname) {
		return(NULL)
})

convertObjectives <- function(flatModel, model) {
	retval <- lapply(flatModel@objectives, function(x) {
		genericObjFunConvert(x, flatModel, model)
	})
	return(retval)
}

translateObjectives <- function(model, namespace) {
	if(is.null(model@objective) && 
		length(omxDependentModels(model)) == 0) {
		return(model)
	}
	flatModel <- omxFlattenModel(model, namespace)
	return(translateObjectivesHelper(model, model, flatModel))
}

translateObjectivesHelper <- function(job, model, flatJob) {
	flatObjective <- flatJob@objectives[[omxIdentifier(model@name, 'objective')]]
	job <- genericObjModelConvert(flatObjective, job, model, flatJob)
	if (length(model@submodels) > 0) {
		for(i in 1:length(model@submodels)) {
			submodel <- model@submodels[[i]]
			if (submodel@independent == FALSE) {
				job <- translateObjectivesHelper(job, submodel, flatJob)
			}
		}
	}
	return(job)
}