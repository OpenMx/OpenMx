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

directories <- c('demo', 'models/passing')

null <- file('/dev/null', 'w')

sink(null, type = 'output')	
sink(null, type = 'message')


files <- list.files(directories, pattern = '.+[.]R',
	full.names = TRUE, recursive = TRUE)
	
errors <- list()

handleError <- function(error) {
	
}
	
errorRecover <- function(script, index) {
	sink(type = 'output')
	sink(type = 'message')	
	cat(paste("Running model", index, "of",
		length(files), script, "...\n"))
	sink(null, type = 'output')	
	sink(null, type = 'message')	
	tryCatch(source(script), error = function(x) {
		errors[[script]] <<- x
	})
}

if (length(files) > 0) {
	for (i in 1:length(files)) {
		errorRecover(files[[i]], i)
	}
}	

sink(type = 'output')
sink(type = 'message')	
close(null)

cat("Number of errors:", length(errors), '\n')
if (length(errors) > 0) {
	fileName <- names(errors)
	for (i in 1:length(errors)) {
		cat("From model", fileName[[i]], ':\n')
		print(errors[[i]])
		cat('\n')
	}
}

cat("Finished testing models.\n")