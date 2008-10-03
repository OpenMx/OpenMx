checkEquals <- function(a, b) {
	if (a != b) {
		stop(paste("Error:", a, "and", b, "are not equal"))
	}	
}

checkTrue <- function(a) {	
	if(!a) {
		stop(paste("Error", match.call()$a, "is not true"))
	}
}