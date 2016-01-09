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

removeTrailingSeparator <- function(x) {
	return(sub('/$', '', x))
}

mxSave <- function(model, chkpt.directory = ".", chkpt.prefix = "") {
	if (!is(model, "MxModel")) {
		stop("'model' argument must be a MxModel object")
	}
	if (!missing(chkpt.directory)) model <- mxOption(model,"Checkpoint Directory", chkpt.directory)
	if (!missing(chkpt.prefix))    model <- mxOption(model,"Checkpoint Prefix", chkpt.prefix)
	model <- mxOption(model,"Checkpoint Units",'evaluations')
	model <- mxOption(model,"Checkpoint Count",1)
	model <- mxModel(model, mxComputeOnce('fitfunction', 'fit'))
	mxRun(model, checkpoint=TRUE, silent=TRUE)
	invisible(TRUE)
}

mxRestore <- function(model, chkpt.directory = ".", chkpt.prefix = "", line=NULL, strict=FALSE) {
	if (!is(model, "MxModel")) {
		stop("'model' argument must be a MxModel object")
	}
	chkpt.directory <- removeTrailingSeparator(chkpt.directory)
	pattern <- paste("^\\Q", chkpt.prefix, "\\E.*(\\.omx)$", sep = '')
	chkpt.files <- list.files(chkpt.directory, full.names = FALSE)
	chkpt.files <- grep(pattern, chkpt.files, perl=TRUE, value=TRUE)
	if(length(chkpt.files) == 0) {
		return(model)
	}
	matchIndex <- match(paste(model$name, 'omx', sep="."), chkpt.files)
	if (strict) {
		if (!is.na(matchIndex)) {
			chkpt.files <- chkpt.files[matchIndex]
		} else {
			stop(paste("Cannot find", omxQuotes(paste(model$name, 'omx', sep=".")),
				   "in", chkpt.directory))
		}
	} else {
		# Move the most likely match to the end so those estimates take precedence.
		chkpt.files <- c(chkpt.files[-matchIndex], paste(model$name, 'omx', sep="."))
	}
	if (length(chkpt.files) > 1 && !is.null(line)) {
		stop(paste("Ambiguous: cannot specify line =", line,
			   "with more than one checkpoint found:",
			   omxQuotes(chkpt.files)))
	}
	if (length(chkpt.files) > 1) {
		message(paste("Loading estimates from more than one checkpoint:",
			      omxQuotes(chkpt.files)))
	}
	allPar <- names(omxGetParameters(model, indep=TRUE, free=NA))
	for(i in 1:length(chkpt.files)) {
		filename <- chkpt.files[[i]]
		filepath <- paste(chkpt.directory, filename, sep = '/')
		checkpoint <- read.table(filepath, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
		ign <- match(c("OpenMxContext","OpenMxNumFree","OpenMxEvals","iterations","timestamp","objective"),
			     colnames(checkpoint))
		ign <- ign[!is.na(ign)]
		toSet <- colnames(checkpoint)[-ign]
		mask <- !is.na(match(toSet, allPar))
		if (all(!mask)) next
		if (is.null(line)) {
			row <- nrow(checkpoint)
		} else {
			row <- line
			if (row < 2 || row > nrow(checkpoint)) {
				warning(paste("Requested line", line,
					      "but checkpoint contains lines 2 to",
					      nrow(checkpoint), "; using the last line"))
				row <- nrow(checkpoint)
			}
		}
		values <- as.numeric(checkpoint[row, -ign])
		model <- omxSetParameters(model, labels=toSet[mask], values=values[mask])
	}
	return(model)
}
