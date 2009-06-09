#
#   Copyright 2007-2009 The OpenMx Project
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


checkEqualDimensions <- function(a, b) {
	if(length(a) == 0 && length(b) == 0) {
		
	} else if (!(length(a) > 0 && length(b) > 0)) {
		stop(paste("One of these is a 0x0 matrix:",
			omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')))) 	
	}
	if((is.vector(a) && length(a) > 1 && !is.vector(b)) || 
		(is.vector(b) && length(b) > 1 && !is.vector(a))) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not both vectors")) 	
	}
	if((is.matrix(a) && (nrow(a) > 1 || ncol(a) > 1) && !is.matrix(b)) || 
		(is.matrix(b) && (nrow(b) > 1 || ncol(b) > 1) && !is.matrix(a))) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not both matrices")) 	
	}	
	if (is.vector(a) && (length(a) != length(b))) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"do not have equal length :",
			length(a), 'and', length(b)))
	}
	if (is.matrix(a) && (nrow(a) > 1 || ncol(a) > 1) && any(dim(a) != dim(b))) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not of equal dimension :", 
			paste(dim(a), collapse = ' x '), 'and', 
			paste(dim(b), collapse = ' x ')))
	}
}

omxCheckEquals <- function(a, b) {
	checkEqualDimensions(a, b)	
	if (any(a != b)) {
		stop(paste(a, "and", b, "are not equal"))
	}	
}

omxCheckSetEquals <- function(a, b) {
	checkEqualDimensions(a, b)	
	if (!setequal(a, b)) {
		stop(paste(a, "and", b, "do not contain the same elements"))
	}	
}

omxCheckTrue <- function(a) {	
	if (any(!a)) {
		stop(paste(match.call()$a, "is not true"))
	}
}


omxCheckCloseEnough <- function(a, b, epsilon=10^(-15)) {
	checkEqualDimensions(a, b)
	check <- any(mapply(function(x,y) {
			abs(x - y) > epsilon }, 
			as.vector(a), as.vector(b)))
	if (check) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"and", omxQuotes(paste(b, collapse = ' ')), 
			"are not equal to within", epsilon))
	}
}

omxCheckWithinPercentError <- function(a, b, epsilon=10^(-15)) {
	checkEqualDimensions(a, b)	
	check <- any(mapply(function(x,y) {
			(abs(x - y)/x * 100) > epsilon }, 
			as.vector(a), as.vector(b)))	
	if (check) {
		stop(paste(omxQuotes(paste(a, collapse = ' ')), 
			"does not estimate", 
			omxQuotes(paste(b, collapse = ' ')), 
			"within", epsilon, "percent"))
	}
}
