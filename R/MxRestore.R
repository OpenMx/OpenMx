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

mxRestore <- function(model, chkpt.directory = ".", chkpt.prefix = "") {	
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
		values <- as.numeric(checkpoint[nrow(checkpoint), -ign])
		model <- omxSetParameters(model, labels=toSet[mask], values=values[mask])
	}
	return(model)
}
