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

removeTrailingSeparator <- function(x) {
	return(sub('/$', '', x))
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
	namespace <- omxGenerateNamespace(model)
	flatModel <- omxFlattenModel(model, namespace)
	pList <- generateParameterList(flatModel)
	for(i in 1:length(chkpt.files)) {
		filename <- chkpt.files[[i]]
		modelname <- substr(filename, nchar(chkpt.prefix) + 1, nchar(filename) - 4)
		filepath <- paste(chkpt.directory, filename, sep = '/')
		checkpoint <- read.table(filepath, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
		model <- restoreCheckpointModel(model, modelname, checkpoint, flatModel, pList)
	}
	return(model)
}

restoreCheckpointModel <- function(model, modelname, checkpoint, flatModel, pList) {
	if (model@independent) {
		namespace <- omxGenerateNamespace(model)
		flatModel <- omxFlattenModel(model, namespace)
		pList <- generateParameterList(flatModel)
	}
	if (modelname == model@name) {
		values <- as.numeric(checkpoint[nrow(checkpoint), 4:ncol(checkpoint)])
		model <- omxUpdateModelValues(model, flatModel, pList, values)
	}
	model@submodels <- lapply(model@submodels, restoreCheckpointModel, modelname, checkpoint, flatModel, pList)
	return(model)
}

