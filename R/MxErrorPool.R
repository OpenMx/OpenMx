#
#   Copyright 2007-2015 The OpenMx Project
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

mxErrorPool <- function(modelnames = NA, reset = FALSE) {
	if (!is.logical(reset) || length(reset) != 1 || is.na(reset)) {
		stop("'reset' argument must be either TRUE or FALSE")
	} else if (!single.na(modelnames) && !is.character(modelnames)) {
		stop("'modelnames' argument must be either NA or a chararacter vector")
	} else if (reset) {
		unlockErrorPool()
		.errorPoolList <<- list()
	}
	if (single.na(modelnames)) {
		retval <- .errorPoolList
	} else {
		retval <- list()
		for(i in 1:length(modelnames)) {
			name <- modelnames[[i]]
			value <- .errorPoolList[[name]]
			if (is.null(value)) {
				retval[name] <- list(NULL)
			} else {
				retval[name] <- value
			}
		}
	}
	return(retval)
}

.errorPoolLocked <- TRUE
.errorPoolList <- list()
