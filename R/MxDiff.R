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

##' Set difference on regular types or S4 objects
##'
##' @param a a
##' @param b b
##' @param slots slots
imxDiff <- function(a, b, slots = c("setequal", "intersect")) {
	if (identical(slots, c("setequal", "intersect"))) {
		slots <- "setequal"
	}
	if (!is.character(slots) || length(slots) != 1 
		|| is.na(slots) || !(slots %in% c("setequal", "intersect"))) {
		stop("'slots' argument must be one of c(\"setequal\", \"intersect\")")
	}
	if (isS4(a) && isS4(b)) {
		diffS4Objects(a, b, slots)
	} else if (is.list(a) && is.list(b)) {
		diffLists(a, b, slots)
	} else {
		stop("arguments 'a' and 'b' must be either both S4 objects or both lists")
	}
}

diffLists <- function(a, b, slots) {
	aSlotNames <- names(a)
	bSlotNames <- names(b)
	if (slots == "setequal" && !setequal(aSlotNames, bSlotNames)) {
		msg <- paste("The 'slots' argument is 'setequal' and the two",
			"lists do not have setequal names.")
		stop(msg)
	}
	if (slots == "setequal") {
		compareNames <- aSlotNames
	} else {
		compareNames <- intersect(aSlotNames, bSlotNames)
	}
	aValues <- a[compareNames]
	bValues <- b[compareNames]
	results <- mapply(identical, aValues, bValues)
	return(compareNames[!results])
}

diffS4Objects <- function(a, b, slots) {
	aSlotNames <- slotNames(class(a))
	bSlotNames <- slotNames(class(b))
	if (slots == "setequal" && !setequal(aSlotNames, bSlotNames)) {
		msg <- paste("The 'slots' argument is 'setequal' and the two",
			"objects do not have setequal slot names.")
		stop(msg)
	}
	if (slots == "setequal") {
		compareNames <- aSlotNames
	} else {
		compareNames <- intersect(aSlotNames, bSlotNames)
	}
	aValues <- lapply(compareNames, function(x) { slot(a, x) })
	bValues <- lapply(compareNames, function(x) { slot(b, x) })
	results <- mapply(identical, aValues, bValues)
	return(compareNames[!results])
}
