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


deparseVector <- function(vec) {
	if (is.character(vec)) {
		inner <- paste(sapply(vec, omxQuotes), collapse = ", ")
	} else {
		inner <- paste(vec, collapse = ", ")
	}
	return(paste("c(", inner, ")", sep = ""))
}

deparseDimnames <- function(dimnames) { 
	if (is.null(dimnames[[1]])) {
		rownames <- "NULL"
	} else {
		rownames <- deparseVector(dimnames[[1]])
	}
	if (is.null(dimnames[[2]])) {
		colnames <- "NULL"
	} else {
		colnames <- deparseVector(dimnames[[2]])
	}
	return(paste("list(", rownames, 
		", ", colnames, ")", sep = ""))
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
		return(paste("matrix(", deparseVector(c(t(object))),
			", nrow = ", nrow(object),
			", ncol = ", ncol(object),
			", byrow = TRUE)", sep = ""))
	}
)

setMethod("mxDeparse", "ZeroMatrix",
	function(object, indent = '   ') {
		return(paste("mxMatrix(type = 'Zero', nrow = ",
			nrow(object), ", ncol = ",
			ncol(object), ", name = ", 
			omxQuotes(object@name), ")", sep = ""))
	}
)

setMethod("mxDeparse", "UnitMatrix",
	function(object, indent = '   ') {
		return(paste("mxMatrix(type = 'Unit', nrow = ",
			nrow(object), ", ncol = ",
			ncol(object), ", name = ", 
			omxQuotes(object@name), ")", sep = ""))
	}
)

setMethod("mxDeparse", "IdenMatrix",
	function(object, indent = '   ') {
		return(paste("mxMatrix(type = 'Iden', nrow = ",
			nrow(object), ", ncol = ",
			ncol(object), ", name = ", 
			omxQuotes(object@name), ")", sep = ""))
	}
)

matrixNAdefault <- function(object, location, retval) {
	if (all(is.na(object))) {		
	} else if (all(apply(object, c(1,2), identical, object[1,1]))) {
		if (is.character(object)) {
			retval <- paste(retval, ", ", location, " = ", 
				omxQuotes(object[1,1]), sep = "")
		} else {
			retval <- paste(retval, ", ", location, " = ", 
				object[1,1], sep = "")
		}
	} else {
		retval <- paste(retval, ", ", location, " = ",
			deparseVector(c(t(object@labels))), sep = "")
	}
	return(retval)
}

setMethod("mxDeparse", "MxMatrix",
	function(object, indent = '   ') {
		type <- sub("Matrix", "", class(object)[[1]], fixed = TRUE)
		retval <- paste("mxMatrix(", omxQuotes(type), sep = "")
		retval <- paste(retval, ", ", nrow(object),
				", ", ncol(object), sep = "")
		if (all(object@free)) {
			retval <- paste(retval, ", free = TRUE", sep = "")
		} else if (any(object@free)) {
			retval <- paste(retval, ", free = ",
				deparseVector(c(t(object@free))), sep = "")
		}
		if (all(object@values == 0)) {
		} else if (all(object@values == object@values[1,1])) {
			retval <- paste(retval, ", values = ", 
				object@values[1,1], sep = "")
		} else {
			retval <- paste(retval, ", values = ",
				deparseVector(c(t(object@values))), sep = "")
		}
		retval <- matrixNAdefault(object@labels, "labels", retval)
		retval <- matrixNAdefault(object@lbound, "lbound", retval)
		retval <- matrixNAdefault(object@ubound, "ubound", retval)
		retval <- paste(retval, ", byrow = TRUE", sep = "")
		if (!is.null(dimnames(object))) {
			retval <- paste(retval, ", dimnames = ", deparseDimnames(dimnames(object)), sep = "")
		}
		retval <- paste(retval, ", name = ", omxQuotes(object@name), sep = "")
		retval <- paste(retval, ')', sep = "")
		return(retval)
	}
)

