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
	chkpt.directory <- removeTrailingSeparator(chkpt.directory)
	chkpt.filename <- paste(chkpt.prefix, model@name, ".omx", sep = '')
	filepath <- paste(chkpt.directory, chkpt.filename, sep = '/')
	print.header <- file.access(filepath) != 0
	pList <- omxGetParameters(model)
	if (length(model@output) == 0) {
		iterations <- 0
		objective <- as.numeric(NA)
	} else {
		iterations <- model@output$iterations
		objective <- model@output$minimum
	}
	timestamp <- date()
	fconnection <- file(filepath, "a")
	if (!isOpen(fconnection, "w")) {
		return(FALSE)
	}
	if (print.header) {
		cat("iterations\t", file=fconnection)
		cat("timestamp\t", file=fconnection)
		cat("objective\t", file=fconnection)
		if (length(pList) > 0) {
			for(i in 1:length(pList)) {
				cat(omxQuotes(names(pList)[[i]]), file=fconnection)
				cat("\t", file=fconnection)
			}
		}
		cat("\n", file=fconnection)
	}
	cat(iterations, file=fconnection)
	cat("\t", file=fconnection)
	cat(omxQuotes(timestamp), file=fconnection)
	cat("\t", file=fconnection)
	cat(objective, file=fconnection)
	cat("\t", file=fconnection)
	if (length(pList) > 0) {
		for(i in 1:length(pList)) {
			cat(pList[[i]], file=fconnection)
			cat("\t", file=fconnection)
		}
	}
	cat("\n", file=fconnection)
	close(fconnection)
	return(TRUE)
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
	namespace <- imxGenerateNamespace(model)
	flatModel <- imxFlattenModel(model, namespace)
	dependencies <- cycleDetection(flatModel)
	dependencies <- transitiveClosure(flatModel, dependencies)
	flatModel <- generateParameterList(flatModel, dependencies, list())
	for(i in 1:length(chkpt.files)) {
		filename <- chkpt.files[[i]]
		modelname <- substr(filename, nchar(chkpt.prefix) + 1, nchar(filename) - 4)
		filepath <- paste(chkpt.directory, filename, sep = '/')
		checkpoint <- read.table(filepath, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
		model <- restoreCheckpointModel(model, modelname, checkpoint, flatModel)
	}
	return(model)
}

restoreCheckpointModel <- function(model, modelname, checkpoint, flatModel) {
	if (model@independent) {
		namespace <- imxGenerateNamespace(model)
		flatModel <- imxFlattenModel(model, namespace)
	}
	if (modelname == model@name) {
		ign <- match(c("OpenMxContext","OpenMxNumFree","OpenMxEvals","iterations","timestamp","objective"),
			     colnames(checkpoint))
		ign <- ign[!is.na(ign)]
		values <- as.numeric(checkpoint[nrow(checkpoint), -ign])
		model <- imxUpdateModelValues(model, flatModel, values)
	}
	model@submodels <- lapply(model@submodels, restoreCheckpointModel, modelname, checkpoint, flatModel)
	return(model)
}

