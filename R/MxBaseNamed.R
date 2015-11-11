#
#   Copyright 2013-2014 The OpenMx Project
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

##' @title MxBaseNamed
##' @name MxBaseNamed-class
##'
##' @description
##' This is an internal class and should not be used directly.  It is the
##' base class for named entities.  Fit functions, expectations, and computes
##' contain this class.
##'
##' @aliases
##' MxBaseNamed
##' @rdname MxBaseNamed-class
setClass(Class = "MxBaseNamed", 
	 representation = representation(
	   name = "character",
	   "VIRTUAL"))

setGeneric("qualifyNames",
	   function(.Object, modelname, namespace) standardGeneric("qualifyNames"))

setGeneric("genericNameToNumber", 
	   function(.Object, flatModel, model) standardGeneric("genericNameToNumber"))

setMethod("genericNameToNumber", signature("MxBaseNamed"),
	  function(.Object, flatModel, model) .Object)

##' imxExtractSlot
##'
##' Checks for and extracts a slot from the object
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param x The object
##' @param name the name of the slot
imxExtractSlot <- function(x, name) {
	if (!.hasSlot(x, name)) {
		return(NULL)
	} else {
		return(slot(x, name))
	}
}

##' imxGetSlotDisplayNames
##'
##' Returns a list of display-friendly object slot names
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param object The object from which to get slot names
##' @param pattern Initial pattern to match (default of '.*' matches any)
##' @param slotList List of slots for which toget display names (default = slotNames(object), i.e., all)
##' @param showDots Include slots whose names start with '.' (default FALSE)
##' @param showEmpty Include slots with length-zero contents (default FALSE)
imxGetSlotDisplayNames <- function(object, pattern='.*', slotList=slotNames(object), showDots=FALSE, showEmpty=FALSE) {
	dotSlots <- slotList[substr(slotList,1,1) == "."]	# Eliminate .<anything> slots
	emptySlots <- slotList[sapply(slotList, 			# Eliminate 0-length slots
					function(x, object) { length(slot(object, x)) == 0 }, 
					object=object)]
	if(!showDots) {
		slotList <- setdiff(slotList, dotSlots)
	}
	if(!showEmpty) {
		slotList <- setdiff(slotList, emptySlots)
	}
	grep(pattern, slotList, value=TRUE)
}

##' imxDefaultGetSlotDisplayNames
##' 
##' Returns a list of display-friendly object slot names
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param x The object from which to get slot names
##' @param pattern Initial pattern to match (default of '.*' matches any)
imxDefaultGetSlotDisplayNames <- function(x, pattern='.*') {
	imxGetSlotDisplayNames(x, pattern)
}


##' imxReplaceSlot
##'
##' Checks for and replaces a slot from the object
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param x object
##' @param name the name of the slot
##' @param value replacement value
##' @param check Check replacement value for validity (default TRUE)
imxReplaceSlot <- function(x, name, value, check=TRUE) {
	if (!.hasSlot(x, name)) {
		stop(paste("Object has no element ", name, ".", sep=""))
        # TODO: Should this return NULL?
	} else {
		slot(x, name, check=check) <- value
		return(x)
	}
}
