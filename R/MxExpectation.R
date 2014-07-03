#
#   Copyright 2007-2014 The OpenMx Project
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
# The virtual base class for all expectations. Expectations contain
# enough information to generate simulated data.
#
setClass(Class = "MxBaseExpectation", 
	 representation = representation(
	   data = "MxCharOrNumber",      # filled in during flattening
	   submodels = "MxOptionalCharOrNumber", # filled in during flattening
	   container = "MxOptionalCharOrNumber", # filled in during flattening
	   "VIRTUAL"),
	 contains = "MxBaseNamed")

setClassUnion("MxExpectation", c("NULL", "MxBaseExpectation"))

setGeneric("genericExpDependencies",
	function(.Object, dependencies) {
	return(standardGeneric("genericExpDependencies"))
})

setGeneric("genericExpRename",
	function(.Object, oldname, newname) {
	return(standardGeneric("genericExpRename"))
})

setGeneric("genericExpFunConvert", 
	function(.Object, flatModel, model, labelsData, defVars, dependencies) {
	return(standardGeneric("genericExpFunConvert"))	
})

setGeneric("genericExpAddEntities",
	function(.Object, job, flatJob, labelsData) {
	return(standardGeneric("genericExpAddEntities"))
})

setGeneric("genericExpConvertEntities",
	function(.Object, flatModel, namespace, labelsData) {
	return(standardGeneric("genericExpConvertEntities"))
})


setMethod("genericExpAddEntities", "MxBaseExpectation",
	function(.Object, job, flatJob, labelsData) {
		return(job)
})

setMethod("genericExpAddEntities", "NULL",
	function(.Object, job, flatJob, labelsData) {
		return(job)
})

setMethod("genericExpConvertEntities", "MxBaseExpectation",
	function(.Object, flatModel, namespace, labelsData) {
		return(flatModel)
})

setMethod("genericExpConvertEntities", "NULL",
	function(.Object, flatModel, namespace, labelsData) {
		return(flatModel)
})

setMethod("genericExpDependencies", "MxBaseExpectation",
	function(.Object, dependencies) {
		return(dependencies)
})

setMethod("genericExpDependencies", "NULL",
	function(.Object, dependencies) {
		return(dependencies)
})

setMethod("genericExpRename", "MxBaseExpectation",
	function(.Object, oldname, newname) {
		return(.Object)
})

setMethod("genericExpRename", "NULL",
	function(.Object, oldname, newname) {
		return(NULL)
})

setMethod("$", "MxBaseExpectation", imxExtractSlot)

setReplaceMethod("$", "MxBaseExpectation",
	function(x, name, value) {
		return(imxReplaceSlot(x, name, value, check=TRUE))
	}
)

setMethod("names", "MxBaseExpectation", slotNames)

convertExpectationFunctions <- function(flatModel, model, labelsData, defVars, dependencies) {
	retval <- lapply(flatModel@expectations, function(ex) {
		ex@container <- imxLocateIndex(flatModel, ex@container, ex@name)
		ex@submodels <- imxLocateIndex(flatModel, ex@submodels, ex@name)
		genericExpFunConvert(ex, flatModel, model, labelsData, defVars, dependencies)
	})
	return(retval)
}

expectationFunctionAddEntities <- function(model, flatModel, labelsData) {

	model@.forcesequential <- FALSE
	model@.newobjects <- FALSE

	expectations <- flatModel@expectations

	if (length(expectations) == 0) {
		return(model)
	}

	for(i in 1:length(expectations)) {
		model <- genericExpAddEntities(expectations[[i]], model, flatModel, labelsData)
	}

	return(model)
}

expectationFunctionConvertEntities <- function(flatModel, namespace, labelsData) {

	expectations <- flatModel@expectations

	if (length(expectations) == 0) {
		return(flatModel)
	}

	for(i in 1:length(expectations)) {
		flatModel <- genericExpConvertEntities(expectations[[i]], flatModel, namespace, labelsData)
	}

	return(flatModel)
}
