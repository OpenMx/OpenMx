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

args <- commandArgs(trailingOnly = TRUE)
if (any(args == 'gctorture')) {
	gctorture(TRUE)
	cat("*** GCTORTURE ENABLED ***\n");
}

library(OpenMx)

options('mxPrintUnitTests' = FALSE)

directories <- c('demo', 'models/passing')
if (any(args == 'nightly')) {
	directories <- c(directories, 'models/nightly')
}

null <- tryCatch(suppressWarnings(file('/dev/null', 'w')),  
	error = function(e) { file('nul', 'w') } )


sink(null, type = 'output')	

if (any(args == 'gctorture')) {
	files <- c('models/passing/rfitfunc.R',
		   'demo/RowObjectiveFIMLBivariateSaturated.R',
		   'models/passing/AlgebraComputePassing.R',
		   'models/passing/TestRowObjective.R',
		   'models/passing/JointFIMLTest.R')
} else {
	files <- list.files(directories, pattern = '^.+[.]R$',
			    full.names = TRUE, recursive = TRUE)
}

if (any(args == 'csolnp')) {
	# remove failing tests using setdiff or whatever TODO
#	files <- c('demo/RowObjectiveSimpleExamples.R')
}

if (any(args == 'lisrel')) {
	files <- grep("LISREL", files, value=TRUE, ignore.case=TRUE)
}

errors <- list()
runtimes <- numeric()

errorRecover <- function(script, index) {
	sink(type = 'output')
	cat(paste("Running model", index, "of",
		length(files), script, "...\n"))
	sink(null, type = 'output')	
	start <- Sys.time()
	tryCatch(source(script, chdir = TRUE), 
		error = function(x) {
			errors[[script]] <<- x
		})
	stop.tm <- Sys.time()
	timeDifference <- stop.tm - start
	runtimes[[script]] <<- as.double(timeDifference, units = "secs")
	rm(envir=globalenv(), 
		list=setdiff(ls(envir=globalenv()), 
			c('errors', 'errorRecover', 'null', 'files', 'directories', 'runtimes')))
}

if (length(files) > 0) {
	for (i in 1:length(files)) {
		errorRecover(files[[i]], i)
	}
}	

sink(type = 'output')
close(null)

cat("Number of errors:", length(errors), '\n')
if (length(errors) > 0) {
	fileName <- names(errors)
	for (i in 1:length(errors)) {
		cat("From model", fileName[[i]], ':\n')
		print(errors[[i]]$message)
		cat('\n')
	}
}

warnings()

write.csv(as.data.frame(runtimes), "runtimes.csv")

cat("Finished testing models.\n")
quit(status=length(errors))
