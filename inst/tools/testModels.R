library(OpenMx)

directories <- c('demo')

files <- list.files(directories, pattern = '.+[.]R',
	full.names = TRUE, recursive = TRUE)
	
errorRecover <- function(script) {
	cat(paste("Running", script, "...\n"))
	sink(file = '/dev/null', type = "output")	
	try(source(script))
	sink()
}
	

invisible(lapply(files, errorRecover))

cat("Finished testing models.\n")