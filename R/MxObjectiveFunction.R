#
#   Copyright 2007-2012 The OpenMx Project
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
        info = "list",
		dependencies = "integer",
		result = "matrix", "VIRTUAL"))

setClassUnion("MxObjective", c("NULL", "MxBaseObjective"))

setGeneric("genericObjFunNamespace", 
	function(.Object, modelname, namespace) {
	return(standardGeneric("genericObjFunNamespace"))
})

setGeneric("genericObjFunConvert", 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
	return(standardGeneric("genericObjFunConvert"))	
})

setGeneric("genericObjAddEntities",
	function(.Object, job, flatJob, labelsData) {
	return(standardGeneric("genericObjAddEntities"))
})

setGeneric("genericObjConvertEntities",
	function(.Object, flatModel, namespace, labelsData) {
	return(standardGeneric("genericObjConvertEntities"))
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

setMethod("genericObjAddEntities", "MxBaseObjective",
	function(.Object, job, flatJob, labelsData) {
		return(job)
})

setMethod("genericObjAddEntities", "NULL",
	function(.Object, job, flatJob, labelsData) {
		return(job)
})

setMethod("genericObjConvertEntities", "MxBaseObjective",
	function(.Object, flatModel, namespace, labelsData) {
		return(flatModel)
})

setMethod("genericObjConvertEntities", "NULL",
	function(.Object, flatModel, namespace, labelsData) {
		return(flatModel)
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

convertObjectiveFunctions <- function(flatModel, model, labelsData, defVars, dependencies) {
	retval <- lapply(flatModel@objectives, genericObjFunConvert, 
		flatModel, model, labelsData, defVars, dependencies)
	return(retval)
}

objectiveFunctionAddEntities <- function(model, flatModel, labelsData) {

	model@.forcesequential <- FALSE
	model@.newobjects <- FALSE

	objectives <- flatModel@objectives

	if (length(objectives) == 0) {
		return(model)
	}

	for(i in 1:length(objectives)) {
		model <- genericObjAddEntities(objectives[[i]], model, flatModel, labelsData)
	}

	return(model)
}


objectiveFunctionModifyEntities <- function(flatModel, namespace, labelsData) {

	objectives <- flatModel@objectives

	if (length(objectives) == 0) {
		return(flatModel)
	}

	for(i in 1:length(objectives)) {
		flatModel <- genericObjConvertEntities(objectives[[i]], flatModel, namespace, labelsData)
	}

	return(flatModel)
}

objectiveReadAttributes <- function(objective, values) {
        attr <- attributes(values)
        attributes(values) <- list('dim' = attr$dim)

		dimnames(values) <- dimnames(objective)
        attr$dim <- NULL

		objective@result <- values
        objective@info <- attr
		return(objective)
}
