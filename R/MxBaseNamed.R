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

setClass(Class = "MxBaseNamed", 
	 representation = representation(
	   name = "character",
	   "VIRTUAL"))

setGeneric("qualifyNames",
	   function(.Object, modelname, namespace) standardGeneric("qualifyNames"))

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
