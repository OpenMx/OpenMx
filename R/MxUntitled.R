
omxUntitledNumber <- function() {
	untitledNumber <<- untitledNumber + 1
	return(untitledNumber)	
}

omxUntitledName <- function() {
	name <- paste("untitled", omxUntitledNumber(), sep="")
	return(name)
}