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
        info = "list",
		result = "matrix", "VIRTUAL"))

setClassUnion("MxObjective", c("NULL", "MxBaseObjective"))

setGeneric("genericObjFunNamespace", 
	function(.Object, modelname, namespace) {
	return(standardGeneric("genericObjFunNamespace"))
})

setGeneric("genericObjFunConvert", 
	function(.Object, flatModel, model, defVars) {
	return(standardGeneric("genericObjFunConvert"))	
})

setGeneric("genericObjModelConvert",
	function(.Object, job, model, namespace, flatJob) {
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
	function(.Object, job, model, namespace, flatJob) {
		job@.newobjects <- FALSE
		job@.newobjective <- FALSE
		job@.newtree <- FALSE
		return(list(job, flatJob))
})

setMethod("genericObjModelConvert", "NULL",
	function(.Object, job, model, namespace, flatJob) {
		job@.newobjects <- FALSE
		job@.newobjective <- FALSE
		job@.newtree <- FALSE
		return(list(job, flatJob))
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

convertObjectives <- function(flatModel, model, defVars) {
	retval <- lapply(flatModel@objectives, function(x) {
		genericObjFunConvert(x, flatModel, model, defVars)
	})
	return(retval)
}

translateObjectives <- function(model, namespace, flatModel) {
	model@.forcesequential <- FALSE
	if(is.null(model@objective) && 
		length(imxDependentModels(model)) == 0) {
		return(list(model, namespace, flatModel))
	}
	return(translateObjectivesHelper(model, namespace, flatModel))
}

translateObjectivesHelper <- function(job, namespace, flatJob) {
	objectives <- flatJob@objectives
	if (length(objectives) == 0) {
		return(list(job, namespace, flatJob))
	}
	for(i in 1:length(objectives)) {
		objective <- objectives[[i]]
		objectivename <- objective@name
		modelname <- unlist(strsplit(objective@name, imxSeparatorChar, fixed=TRUE))[[1]]
		model <- job[[modelname]]
		repeat {
			if(is.null(objective)) { break }
			# .newobjects is TRUE when new matrices or algebras are created
			# .newobjective is TRUE when one objective function is transformed into another
			# .newtree is TRUE when one objective function is transformed into multiple objective functions
			job@.newobjects <- TRUE
			job@.newobjective <- TRUE
			job@.newtree <- TRUE
			pair <- genericObjModelConvert(objective, job, model, namespace, flatJob)
			job <- pair[[1]]
			flatJob <- pair[[2]]
			if (job@.newtree) {
				namespace <- imxGenerateNamespace(job)
				flatJob <- imxFlattenModel(job, namespace)
				return(translateObjectivesHelper(job, namespace, flatJob))
			}
			if (job@.newobjects) {
				namespace <- imxGenerateNamespace(job)
				flatJob <- imxFlattenModel(job, namespace)
			}
			if (job@.newobjective) {
				model <- job[[modelname]]
				objective <- flatJob[[objectivename]]
			} else { break }
		}
	}
	return(list(job, namespace, flatJob))
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
