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

omxLapply <- function(x, fun, ...) {
	if (length(x) == 0) return(x)
	libraries <- search()
	if ("package:snowfall" %in% libraries) {
		return(sfLapply(x, fun, ...))
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
