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

##' imxMpiWrap
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param fun fun
imxMpiWrap <- function(fun) {
	function(...) {
		require(OpenMx, quietly = TRUE)
		return(fun(...))
	}
}

##' imxSfClient
##'
##' As of snowfall 1.84, the snowfall supervisor process
##' stores an internal state information in a variable 
##' named ".sfOption" that is located in the "snowfall" namespace.
##' The snowfall client processes store internal state
##' information in a variable named ".sfOption" that is located
##' in the global namespace.
##' 
##' As long as the previous statement is true, then the current
##' process is a snowfall client if-and-only-if exists(".sfOption").
imxSfClient <- function() {
	return(exists(".sfOption"))
}

omxLapply <- function(x, fun, ...) {
	if (length(x) == 0) return(x)
	libraries <- search()
	if ("package:snowfall" %in% libraries) {
		return(snowfall::sfClusterApplyLB(x, fun, ...))
	} else if ("package:Rmpi" %in% libraries) {
		return(Rmpi::mpi.parLapply(x, imxMpiWrap(fun), ...))
	} else {
		return(lapply(x, fun, ...))
	}
}

omxSapply <- function(x, fun, ..., simplify = TRUE, USE.NAMES = TRUE) {
	if (length(x) == 0) return(x)
	libraries <- search()
	if ("package:snowfall" %in% libraries) {
		return(snowfall::sfSapply(x, fun, ..., simplify = simplify, USE.NAMES = USE.NAMES))
	} else {
		return(sapply(x, fun, ..., simplify = simplify, USE.NAMES = USE.NAMES))
	}
}

omxApply <- function(x, margin, fun, ...) {
	if (length(x) == 0) return(x)
	libraries <- search()
	if ("package:snowfall" %in% libraries) {
		return(snowfall::sfApply(x, margin, fun, ...))
	} else {
		return(apply(x, margin, fun, ...))
	}
}
