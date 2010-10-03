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

unlockErrorPool <- function() {
	if(.errorPoolLocked) {
		namespace <- getNamespace("OpenMx")
		unlockBinding('.errorPoolLocked', namespace)
		.errorPoolLocked <<- FALSE
		lockBinding('.errorPoolLocked', namespace)
		unlockBinding('.errorPoolList', namespace)
	}
}

mxErrorPool <- function(..., reset = FALSE) {
	if (!is.logical(reset) || length(reset) != 1 || is.na(reset)) {
		stop("'reset' argument must be either TRUE or FALSE")
	} else if (reset) {
		unlockErrorPool()
		.errorPoolList <<- list()
	}
	args <- list(...)
	if (length(args) == 0) {
		retval <- .errorPoolList	
	} else {
		allchars <- sapply(args, is.character)
		if (!all(allchars)) {
			stop("All '...' values to mxErrorPool must be character vectors")
		}
		args <- unlist(args, use.names = FALSE)
		retval <- list()
		for(i in 1:length(args)) {
			arg <- args[[i]]
			value <- .errorPoolList[[arg]]
			if (is.null(value)) {
				retval[arg] <- list(NULL)
			} else {
				retval[arg] <- value
			}
		}
	}
	return(retval)
}

.errorPoolLocked <- TRUE
.errorPoolList <- list()
