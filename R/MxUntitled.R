
omxUntitledNumber <- function() {
	omxUntitledCounter <<- omxUntitledCounter + 1
	return(omxUntitledCounter)	
}

omxUntitledName <- function() {
	name <- paste("untitled", omxUntitledNumber(), sep="")
	return(name)
}

omxUntitledCounter <- 0
