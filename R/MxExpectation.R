#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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
##' MxBaseExpectation
##'
##' @description
##' The virtual base class for all expectations. Expectations contain
##' enough information to generate simulated data.  This is an internal
##' class and should not be used directly.
##'
##' @aliases
##' $,MxBaseExpectation-method
##' $<-,MxBaseExpectation-method
##' @seealso
##' \link{mxExpectationNormal}, \link{mxExpectationRAM}, \link{mxExpectationLISREL}, \link{mxExpectationStateSpace},
##' \link{mxExpectationBA81}
##' @rdname MxBaseExpectation-class
setClass(Class = "MxBaseExpectation", 
	 representation = representation(
		 data = "MxCharOrNumber",      # filled in during flattening
		 dataColumns = "integer",      # deprecated
		 dataColumnNames = "MxOptionalChar",  # subset and permutation of data columns
	     .runDims = "character",
	     output = "list",
	     debug = "list",
	   "VIRTUAL"),
	 contains = "MxBaseNamed")

##' @title MxExpectation
##' @name MxExpectation-class
##'
##' @description
##' This is an internal class and should not be used directly.
##'
##' @aliases
##' MxExpectation
##' MxExpectation-class
##' @rdname MxExpectationXclass
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
	function(.Object, flatModel, model, labelsData, dependencies) {
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

setGeneric("genericGetExpected",
	function(.Object, model, what, defvar.row, subname) {
	return(standardGeneric("genericGetExpected"))
})

setGeneric("genericGetExpectedVector",
	function(.Object, model, defvar.row, subname) {
	return(standardGeneric("genericGetExpectedVector"))
})

setGeneric("genericGetExpectedStandVector",
	function(.Object, model, defvar.row, subname) {
	return(standardGeneric("genericGetExpectedStandVector"))
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
		stop("Not implemented yet")
})

setMethod("genericExpRename", "NULL",
	function(.Object, oldname, newname) {
		return(NULL)
})

setMethod("genericGetExpected", "MxBaseExpectation",
	function(.Object, model, what, defvar.row, subname) stop("Not implemented"))

setMethod("genericGetExpectedVector", "MxBaseExpectation",
	function(.Object, model, defvar.row, subname) stop("Not implemented"))

setMethod("genericGetExpectedStandVector", "MxBaseExpectation",
	function(.Object, model, defvar.row, subname) stop("Not implemented"))

setMethod("$", "MxBaseExpectation", imxExtractSlot)

setReplaceMethod("$", "MxBaseExpectation",
	function(x, name, value) {
		return(imxReplaceSlot(x, name, value, check=TRUE))
	}
)

setMethod("names", "MxBaseExpectation", slotNames)

convertExpectationFunctions <- function(flatModel, model, labelsData, dependencies) {
	# The idea is to split genericExpFunConvert into sanity checking,
	# which often requires access to the actual objects, and
	# converting symbolic names into numbers that are easy to deal
	# with in the backend.
	retval <- lapply(flatModel@expectations, function(ex) {
		genericExpFunConvert(ex, flatModel, model, labelsData, dependencies)
	})
	retval <- lapply(retval, function(ex) {
		genericNameToNumber(ex, flatModel, model)
	})
	return(retval)
}

expectationFunctionAddEntities <- function(model, flatModel, labelsData) {

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

#By the time this function is called, the optionsList will contain any locally set mxOptions:
getPrecisionPerExpectation <- function(expectation, optionsList){
	
	#Which options need a proper numerical value?
	needStepSize <- as.logical(length(grep(pattern="Auto",x=optionsList$"Gradient step size",ignore.case=T)))
	needIters <- as.logical(length(grep(pattern="Auto",x=optionsList$"Gradient iterations",ignore.case=T)))
	needFuncPrec <- as.logical(length(grep(pattern="Auto",x=optionsList$"Function precision",ignore.case=T)))
	
	#Do we have an ordinal-thresholds expectation?:
	isOrdThresh <- (class(expectation) %in% c("MxExpectationNormal","MxExpectationLISREL","MxExpectationRAM")) && 
		!single.na(expectation@thresholds)
	
	if(needStepSize){
		if(isOrdThresh){stepSize <- 1.0e-5} #<--Default value, times 1e3
		else{stepSize <- 1.0e-7} #<--Default value
	} 
	else{stepSize <- as.numeric(optionsList$"Gradient step size")}
	
	if(needIters){
		#Neither 3L nor 4L agrees with the default value of 1L...but this how the definitions for 
		#generic method "genericExpGetPrecision" were written...:
		if(isOrdThresh){iterations <- 3L}
		else{iterations <- 4L}
	}
	else{iterations <- as.integer(optionsList$"Gradient iterations")}
	
	if(needFuncPrec){
		if(isOrdThresh){functionPrecision <- 1e-10}
		else{functionPrecision <- 1e-14} #<--Default value
	}
	else{functionPrecision <- as.numeric(optionsList$"Function precision")}
	
	return(list(stepSize=stepSize, iterations=iterations, functionPrecision=functionPrecision))
}


##' imxHasThresholds
##'
##' This is an internal function exported for those people who know
##' what they are doing.  This function checks if a model (or its
##' submodels) has any thresholds.
##'
##' @param model model
imxHasThresholds <- function(model) {
	if(length(model@expectation) && 
		 (class(model@expectation) %in% c("MxExpectationNormal","MxExpectationLISREL","MxExpectationRAM")) && 
		 !single.na(model@expectation@thresholds)  ){
		return(TRUE)
	}
	# Check submodels
	if(length(model@submodels)){
		attempt <- sapply(model@submodels, imxHasThresholds)
		if(any(attempt)){
			return(TRUE)
		}
	}
	return(FALSE)
}
