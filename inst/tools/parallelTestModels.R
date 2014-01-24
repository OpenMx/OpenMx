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

library(OpenMx)
library(snowfall)

sfInit(parallel=TRUE)

sfLibrary(snowfall)
sfLibrary(OpenMx)

options('mxPrintUnitTests' = FALSE)

directories <- c('demo', 'models/passing')

files <- list.files(directories, pattern = '^.+[.]R$',
	full.names = TRUE, recursive = TRUE)
	
errorRecover <- function(index) {
	retval <- list()
	script <- files[[index]]
	cat(paste("Running model", index, "of",
		length(files), script, "...\n"))
	tryCatch(source(script, chdir = TRUE), 
		error = function(x) {
			retval[[script]] <<- x
		})
	rm(envir=globalenv(),
                list=setdiff(ls(envir=globalenv()),
                        c('errorRecover', 'files', 'directories')))
	return(retval)
}

sfExportAll()

if (length(files) > 0) {
	results <- sfLapply(1:length(files), errorRecover)
	errors <- unlist(results, recursive=FALSE)
}	

cat("Number of errors:", length(errors), '\n')
if (length(errors) > 0) {
	fileName <- names(errors)
	for (i in 1:length(errors)) {
		cat("From model", fileName[[i]], ':\n')
		print(errors[[i]]$message)
		cat('\n')
	}
}

cat("Finished testing models.\n")

sfStop()
