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

generateCommunicationList <- function(model, checkpoint, useSocket, options) {
	if (mxOption(model,'Always Checkpoint') == "Yes") {
		checkpoint <- TRUE
	}
	if (!checkpoint) return(list())
	retval <- list()
	if (checkpoint) {		
		chkpt.directory <- mxOption(model, 'Checkpoint Directory')
		chkpt.directory <- removeTrailingSeparator(chkpt.directory)
		chkpt.prefix <- mxOption(model, 'Checkpoint Prefix')
		chkpt.units <- mxOption(model, 'Checkpoint Units')
		chkpt.count <- mxOption(model, 'Checkpoint Count')
		chkpt.count <- mxOption(model, 'Checkpoint Count')
		if (length(chkpt.count) == 2) {
			chkpt.count <- chkpt.count[[chkpt.units]]
		}
		if (is.null(chkpt.count)) chkpt.count <- .Machine$integer.max

		if (!is.numeric(chkpt.count) || chkpt.count < 0) {
			stop(paste("'Checkpoint Count' model option",
				"must be a non-negative value in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (!(is.character(chkpt.prefix) && length(chkpt.prefix) == 1)) {
			stop(paste("'Checkpoint Prefix' model option",
				"must be a string in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (!(is.character(chkpt.directory) && length(chkpt.directory) == 1)) {
			stop(paste("'Checkpoint Directory' model option",
				"must be a string in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (!(is.character(chkpt.units) && length(chkpt.units) == 1)) {
			stop(paste("'Checkpoint Units' model option",
				"must be a string in", 
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		filename <- paste(chkpt.prefix, paste(model$name, 'omx', sep = '.'), sep = '')
		fullpath <- paste(chkpt.directory, filename, sep="/")
		override <- mxOption(model, "Checkpoint Fullpath")
		if (nchar(override)) {
			fullpath <- override
		}
		description <- list(0L, fullpath, chkpt.units, chkpt.count)
		retval[[length(retval) + 1]] <- description
	}
	return(retval)
}
