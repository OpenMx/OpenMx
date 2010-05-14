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

generateCommunicationList <- function(modelname, chkpt.directory, chkpt.prefix,
		chkpt.units, chkpt.count, sock.server, sock.port, sock.units, sock.count) {
	retval <- list()
	chkpt.directory <- removeTrailingSeparator(chkpt.directory)
	if (!is.numeric(chkpt.count) || chkpt.count < 0) {
		stop(paste("'chkpt.count' argument to mxRun",
			"must be a non-negative value in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if (!is.numeric(sock.count) || sock.count < 0) {
		stop(paste("'sock.count' argument to mxRun",
			"must be a non-negative value in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if ((!is.null(sock.server) && is.null(sock.port)) ||
		(is.null(sock.server) && !is.null(sock.port))) {
		stop(paste("Both 'sock.server' and 'sock.port'",
			"must be specified in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if ((chkpt.count == 0) && (!identical(chkpt.directory, ".") ||
		!identical(chkpt.prefix, "") || !identical(chkpt.units, c("minutes", "iterations")))) {
		stop(paste("'chkpt.count' argument to mxRun",
			"must be specified to use checkpoints in", 
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
	if (chkpt.count > 0) {
		if (!(is.character(chkpt.prefix) && length(chkpt.prefix) == 1)) {
			stop(paste("'chkpt.prefix' argument to mxRun",
				"must be a string in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (!(is.character(chkpt.directory) && length(chkpt.directory) == 1)) {
			stop(paste("'chkpt.directory' argument to mxRun",
				"must be a string in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (identical(chkpt.units, c("minutes", "iterations"))) {
			chkpt.units <- "minutes"
		}
		if (!(is.character(chkpt.units) && length(chkpt.units) == 1)) {
			stop(paste("'chkpt.units' argument to mxRun",
				"must be a string in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		filename <- paste(chkpt.prefix, paste(modelname, 'omx', sep = '.'), sep = '')
		if (chkpt.units == "minutes") {
			type <- 0L
		} else if (chkpt.units == "iterations") {
			type <- 1L
		} else {
			stop(paste("'chkpt.units' argument to mxRun",
				"must be either 'minutes' or 'iterations' in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		description <- list(0L, chkpt.directory, filename, type, chkpt.count)
		retval[[length(retval) + 1]] <- description
	}
	if (!is.null(sock.server) && !is.null(sock.port)) {
		if (!(is.character(sock.server) && length(sock.server) == 1)) {
			stop(paste("'sock.server' argument to mxRun",
				"must be a string in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (!(is.numeric(sock.port) && length(sock.port) == 1)) {
			stop(paste("'sock.port' argument to mxRun",
				"must be a numeric value in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (identical(sock.units, c("minutes", "iterations"))) {
			sock.units <- "minutes"
		}
		if (!(is.character(sock.units) && length(sock.units) == 1)) {
			stop(paste("'sock.units' argument to mxRun",
				"must be a string in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (sock.count == 0) {
			stop(paste("'sock.count' argument to mxRun",
				"must be a non-negative value in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (sock.units == "minutes") {
			type <- 0L
		} else if (sock.units == "iterations") {
			type <- 1L
		} else {
			stop(paste("'sock.units' argument to mxRun",
				"must be either 'minutes' or 'iterations' in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		description <- list(1L, sock.server, as.integer(sock.port), type, sock.count)
		retval[[length(retval) + 1]] <- description
	}
	return(retval)
}
