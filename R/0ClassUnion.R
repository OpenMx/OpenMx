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


##' A character or integer
##' @name MxCharOrNumber-class
setClassUnion("MxCharOrNumber", c("character", "integer"))

##' An optional character
##' @name MxOptionalChar-class
setClassUnion("MxOptionalChar", c("NULL", "character"))

##' An optional data.frame
##' @name MxOptionalDataFrame-class
setClassUnion("MxOptionalDataFrame", c("NULL", "data.frame"))

##' @title An optional logical
##' @name MxOptionalLogical-class
##' @rdname MxOptionalLogical-class
##' @description This is an internal class, the union of NULL and logical.
setClassUnion("MxOptionalLogical", c("NULL", "logical"))

##' A character, integer, or NULL
##' @name MxOptionalCharOrNumber-class
setClassUnion("MxOptionalCharOrNumber", c("NULL", "character", "integer"))

##' An optional list
##' @name MxListOrNull-class
setClassUnion("MxListOrNull", c("list", "NULL"))

##' A character, list or NULL
##' @name MxCharOrList-class
setClassUnion("MxCharOrList", c("character", "list"))

##' An optional matrix
##' @name MxOptionalMatrix-class
setClassUnion("MxOptionalMatrix", c("NULL", "matrix"))

##' An optional numeric
##' @name MxOptionalNumeric-class
setClassUnion("MxOptionalNumeric", c("NULL", "numeric"))

##' An optional integer
##' @name MxOptionalInteger-class
setClassUnion("MxOptionalInteger", c("NULL", "integer"))

##' A character or logical
##' @name MxCharOrLogical-class
setClassUnion("MxCharOrLogical", c("character", "logical"))

##' A package_version or character
##' @name MxVersionType-class
setOldClass('package_version')
setClassUnion("MxVersionType", c("package_version", "character"))

factorize <- function(x, levels, labels, exclude, collapse) {
	x <- as.character(x)
	if (length(exclude) && all(!is.na(exclude))) {
		overlap <- match(exclude, levels)
		if (any(!is.na(overlap))) {
			msg <- paste("Factor levels and exclude vector are not disjoint; both contain",
				     omxQuotes(levels[overlap]))
			stop(msg)
		}
		x[which(!is.na(match(x, exclude)))] <- NA
	}
	noMatch <- !is.na(x) & is.na(match(x, levels))
	if (any(noMatch)) {
		msg <- paste("The following values are not mapped to factor levels and not excluded:",
			     omxQuotes(unique(x[noMatch])))
		stop(msg)
	}
	if (collapse) {
    corder <- order(labels)
    cLabels <- labels[corder]
    cLevels <- levels[corder]
	  dups <- duplicated(cLabels)
	  newLevels <- cLevels[!dups]
		notDup <- which(!dups)
		for (dx in which(dups)) {
			from <- cLevels[dx]
			to <- newLevels[findInterval(dx, notDup)]
			x[x==from] <- to
		}
    mask <- !duplicated(labels)
		levels <- levels[mask]
		labels <- labels[mask]
	} else {
	  dups <- duplicated(labels)
	  if (any(dups)) stop(paste("Duplicate labels and collapse=TRUE not specified:",
					  omxQuotes(unique(labels[dups]))))
	}

	f <- factor(x, levels, labels, exclude, ordered=TRUE)
	attr(f, 'mxFactor') <- TRUE
	f
}

mxFactor <- function(x = character(), levels, labels = levels, exclude = NA, ordered = TRUE, collapse=FALSE) {
	if(missing(levels)) {
		stop("the 'levels' argument is not optional")
	}
	if(!identical(ordered, TRUE)) {
		stop("the 'ordered' argument must be TRUE")
	}
	if (is.data.frame(x)) {
		if (is.list(levels)) {
			return(data.frame(mapply(factorize, x, levels, labels,
				MoreArgs=list(exclude = exclude, collapse=collapse), SIMPLIFY=FALSE),
				check.names = FALSE, row.names=rownames(x)))
		} else {
			return(data.frame(lapply(x, factorize, levels, labels, exclude, collapse),
				check.names = FALSE, row.names=rownames(x)))
		}
	} else if (is.matrix(x)) {
		stop(paste("argument 'x' to mxFactor()",
		"is of illegal type matrix,",
		"legal types are vectors or data.frames"))
	} else {
		return(factorize(x, levels, labels, exclude, collapse))
	}
}

prohibitDotdotdot <- function(args) {
  if (length(args) == 0) return()
  stop(paste0(as.character(sys.call(-1))[1], " does not accept ... arguments. ",
              "The first parameter in ... was named ", omxQuotes(names(args)[1]), 
              " with value '", args[[1]], "'"),
       call.=FALSE)
}

