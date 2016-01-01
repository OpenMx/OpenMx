#
#   Copyright 2007-2016 The OpenMx Project
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

componentCombination <- function(func, slotname, args) {
	args <- lapply(args, slot, slotname)
	return(do.call(func, args))
}

checkAnonymousFreeParameters <- function(mxmatrix) {
	if(any(is.na(mxmatrix@labels) & mxmatrix@free)) {
		msg <- paste("Anonymous free parameters are",
			"not allowed when 'allowUnlabeled' argument",
			"is FALSE.")
		stop(msg, call. = FALSE)
	}
}

omxCbind <- function(..., allowUnlabeled = 
		getOption("mxOptions")[["Allow Unlabeled"]], 
		dimnames = NA, name = NA) {

	args <- list(...)
	objcheck <- sapply(args, isS4)
	if (any(!objcheck)) {
		stop("All input arguments must be MxMatrix objects.")
	}
	objcheck <- sapply(args, is, "MxMatrix")
	if (any(!objcheck)) {
		stop("All input arguments must be MxMatrix objects.")
	}
	if (length(allowUnlabeled) != 1 || 
		!is.logical(allowUnlabeled) || 
		is.na(allowUnlabeled)) {
		stop("'allowUnlabeled' must be either TRUE or FALSE.")
	}
	if (!allowUnlabeled) {
		lapply(args, checkAnonymousFreeParameters)
	}
	values <- componentCombination(cbind, "values", args)
	free   <- componentCombination(cbind, "free", args)
	labels <- componentCombination(cbind, "labels", args)
	lbound <- componentCombination(cbind, "lbound", args)
	ubound <- componentCombination(cbind, "ubound", args)
	retval <- mxMatrix(type = "Full", free = free, 
		values = values, labels = labels, lbound = lbound,
		ubound = ubound, dimnames = dimnames, name = name)
	return(retval)
}

omxRbind <- function(..., allowUnlabeled = 
		getOption("mxOptions")[["Allow Unlabeled"]], 
		dimnames = NA, name = NA) {

	args <- list(...)
	objcheck <- sapply(args, isS4)
	if (any(!objcheck)) {
		stop("All input arguments must be MxMatrix objects.")
	}
	objcheck <- sapply(args, is, "MxMatrix")
	if (any(!objcheck)) {
		stop("All input arguments must be MxMatrix objects.")
	}
	if (length(allowUnlabeled) != 1 || 
		!is.logical(allowUnlabeled) || 
		is.na(allowUnlabeled)) {
		stop("'allowUnlabeled' must be either TRUE or FALSE.")
	}
	if (!allowUnlabeled) {
		lapply(args, checkAnonymousFreeParameters)
	}
	values <- componentCombination(rbind, "values", args)
	free   <- componentCombination(rbind, "free", args)
	labels <- componentCombination(rbind, "labels", args)
	lbound <- componentCombination(rbind, "lbound", args)
	ubound <- componentCombination(rbind, "ubound", args)
	retval <- mxMatrix(type = "Full", free = free, 
		values = values, labels = labels, lbound = lbound,
		ubound = ubound, dimnames = dimnames, name = name)
	return(retval)
}

omxTranspose <- function(matrix, allowUnlabeled = 
		getOption("mxOptions")[["Allow Unlabeled"]], 
		dimnames = NA, name = NA) {

	if (!isS4(matrix)) {
		stop("input argument must be a MxMatrix object.")
	}
	if (!is(matrix, "MxMatrix")) {
		stop("input argument must be a MxMatrix object.")
	}
	if (length(allowUnlabeled) != 1 || 
		!is.logical(allowUnlabeled) || 
		is.na(allowUnlabeled)) {
		stop("'allowUnlabeled' must be either TRUE or FALSE.")
	}
	if (!allowUnlabeled) {
		checkAnonymousFreeParameters(matrix)
	}
	values <- t(matrix@values)
	free   <- t(matrix@free)
	labels <- t(matrix@labels)
	lbound <- t(matrix@lbound)
	ubound <- t(matrix@ubound)
	retval <- mxMatrix(type = "Full", free = free, 
		values = values, labels = labels, lbound = lbound,
		ubound = ubound, dimnames = dimnames, name = name)
}
