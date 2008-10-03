setClass(Class = "MxSparseMatrix",
	representation = representation(
		rowVector = "vector",
		colVector = "vector",
		dataVector = "vector",
		nrow = "numeric",
		ncol = "numeric"),
	prototype = prototype(
		rowVector = vector(mode="numeric",0),
		colVector = vector(mode="numeric",0),
		dataVector = vector(mode="character",0),
		nrow = 1, ncol = 1))

setMethod("initialize", "MxSparseMatrix",
	function(.Object, data = NA, nrow = 1, ncol = 1) {
		if (is.matrix(data)) {
			tdata <- t(data)
			nonZero <- (data != 0) | is.na(data)
			rowMatrix <- t(row(data) * nonZero)
			colMatrix <- t(col(data) * nonZero)
			rowVector <- as.vector(rowMatrix)
			colVector <- as.vector(colMatrix)
			.Object@dataVector <- tdata[tdata != 0]
			.Object@rowVector <- rowVector[rowVector > 0]
			.Object@colVector <- colVector[colVector > 0]
			.Object@nrow <- nrow(data)
			.Object@ncol <- ncol(data)
		} else {
			.Object@nrow <- nrow
			.Object@ncol <- ncol
	    }
		return(.Object)
	}
)					

setMethod("show", "MxSparseMatrix",
	function(object) {
	    cat(paste(nrow(object), "x", ncol(object), "MxSparseMatrix matrix \n"))
		print(as.matrix(object))
	}
)

setMethod("t", "MxSparseMatrix",		
	function(x) {
	    return(new("MxSparseMatrix", data=t(as.matrix(x))))
		return(result)
	}
)

setMethod("nnzero", "MxSparseMatrix",
	function(x, na.counted = NA) {
		return(length(x@rowVector))
	}
)
		
setMethod("as.matrix", "MxSparseMatrix",
	function(x, ...) {
		result <- matrix(0, nrow(x), ncol(x))
		result[x@rowVector + (x@colVector - 1)*nrow(x)] <- x@dataVector
		return(result)
	}
)

setMethod("nrow", "MxSparseMatrix",
	function(x) {
	    return(x@nrow)
	}
)

setMethod("ncol", "MxSparseMatrix",
	function(x) {
	    return(x@ncol)
	}
)



setMethod("[", "MxSparseMatrix",
	function(x, i, j, ..., drop = FALSE) {
	    if (i > x@nrow || j > x@ncol) {
	        stop("subscript out of bounds")
	    } else if (length(x@rowVector) == 0) {
			return(0)
        } else {
			ubound <- findInterval(i, x@rowVector)
			if ((ubound > 0) && x@rowVector[ubound] == i) {
				lbound <- findInterval(i - 1, x@rowVector) + 1
				index <- findInterval(j, x@colVector[lbound:ubound])
				offset <- index + lbound - 1
				if (offset > 0 && x@colVector[offset] == j) {
					return(x@dataVector[offset])
				}
			} 
		}
		return(0)
    }
)

setReplaceMethod("[", "MxSparseMatrix", 
	function(x, i, j, value) {
	    if (i > x@nrow || j > x@ncol) {
	        stop("subscript out of bounds")
        } else if (length(x@rowVector) == 0) {
			x@rowVector <- c(i)
			x@colVector <- c(j)
			x@dataVector <- c(value)
        } else {
			ubound <- findInterval(i, x@rowVector)
			if ((ubound > 0) && x@rowVector[ubound] == i) {
				lbound <- findInterval(i - 1, x@rowVector) + 1
				index <- findInterval(j, x@colVector[lbound:ubound])
				offset <- index + lbound - 1
				if (offset > 0 && x@colVector[offset] == j) {
				    if (value == 0) {
					   x@rowVector <- x@rowVector[-offset]
				       x@colVector <- x@colVector[-offset]
				       x@dataVector <- x@dataVector[-offset]
				    } else {
  					   x@rowVector[offset] <- i
					   x@colVector[offset] <- j
					   x@dataVector[offset] <- value
				    }
				} else if (value != 0) {
					x@rowVector <- append(x@rowVector, i, after = offset)
					x@colVector <- append(x@colVector, j, after = offset)
					x@dataVector <- append(x@dataVector, value, after = offset)
				}
			} else if (value != 0) {
				x@rowVector <- append(x@rowVector, i, after = ubound)
				x@colVector <- append(x@colVector, j, after = ubound)
				x@dataVector <- append(x@dataVector, value, after = ubound)
			}
		}
		return(x)
    }
)
