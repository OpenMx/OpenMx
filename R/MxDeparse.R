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


omxDeparseVector <- function(vec) {
	return(paste("c(",
		paste(vec, collapse = ", "),
		")", sep = ""))
}

setGeneric("mxDeparse", function(object, indent = '   ') { 
	return(standardGeneric("mxDeparse")) 
})

setMethod("mxDeparse", "MxAlgebra",
	function(object, indent = '   ') {
		return(paste("mxAlgebra", "(", 
			deparse(object@formula), ", name = ",
			omxQuotes(object@name), ")", sep = ""))
	}
)

setMethod("mxDeparse", "MxConstraint",
	function(object, indent = '   ') {
		return(paste("mxConstraint", "(",
			omxQuotes(object@alg1), ", ",
			omxQuotes(object@relation), ", ",
			omxQuotes(object@alg2), ", ",
			"name=", omxQuotes(object@name), 
			")", sep = ""))
	}
)

setMethod("mxDeparse", "MxData",
	function(object, indent = '   ') {
		if (is.na(object@means)) {
			means <- "NA"
		} else {
			means <- "..."
		}
		if (is.na(object@numObs)) {
			numObs <- as.character(object@numObs)
		} else {
			numObs <- "NA"
		}
		return(paste("mxData", "(",
			"observed = ...", ", ",
			"type = ", omxQuotes(object@type), ", ",
			"means = ", means, ", ",
			"numObs = ", numObs, ")", sep = ""))
	}
)

setMethod("mxDeparse", "matrix",
	function(object, indent = '   ') {
		if (nrow(object) == 0 || ncol(object) == 0) {
			return(paste("matrix(nrow = ",
				nrow(object), ", ncol = ",
				ncol(object), ")", sep = ""))
		}
		first <- object[1,1]
		if (all(apply(object, c(1,2), identical, first))) {
			return(paste("matrix(", first,
				", nrow = ", nrow(object),
				", ncol = ", ncol(object), ")", sep = ""))
		}
		return(paste("matrix(", omxDeparseVector(c(t(object))),
			", nrow = ", nrow(object),
			", ncol = ", ncol(object),
			", byrow = TRUE)", sep = ""))
	}
)