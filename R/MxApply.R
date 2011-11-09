#
#   Copyright 2007-2011 The OpenMx Project
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

imxMpiWrap <- function(fun) {
	function(...) {
		require(OpenMx, quietly = TRUE)
		return(fun(...))
	}
}

imxSfClient <- function() {
	return("snowfall" %in% loadedNamespaces())
}

omxLapply <- function(x, fun, ...) {
	if (length(x) == 0) return(x)
	libraries <- search()
	if ("package:Swift" %in% libraries && identical(fun,mxRun)) {
		return(swiftLapply(x, fun, ...))
	} else if ("package:snowfall" %in% libraries) {
		return(sfClusterApplyLB(x, fun, ...))
	} else if ("package:Rmpi" %in% libraries) {
		return(mpi.parLapply(x, imxMpiWrap(fun), ...))
	} else {
		return(lapply(x, fun, ...))
	}
}

omxSapply <- function(x, fun, ..., simplify = TRUE, USE.NAMES = TRUE) {
	if (length(x) == 0) return(x)
	libraries <- search()
	if ("package:snowfall" %in% libraries) {
		return(sfSapply(x, fun, ..., simplify = simplify, USE.NAMES = USE.NAMES))
	} else {
		return(sapply(x, fun, ..., simplify = simplify, USE.NAMES = USE.NAMES))
	}
}

omxApply <- function(x, margin, fun, ...) {
	if (length(x) == 0) return(x)
	libraries <- search()
	if ("package:snowfall" %in% libraries) {
		return(sfApply(x, margin, fun, ...))
	} else {
		return(apply(x, margin, fun, ...))
	}
}
