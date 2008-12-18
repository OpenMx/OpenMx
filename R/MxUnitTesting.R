omxCheckEquals <- function(a, b) {
	if (a != b) {
		stop(paste("Error:", a, "and", b, "are not equal"))
	}	
}

omxCheckTrue <- function(a) {	
	if(!a) {
		stop(paste("Error", match.call()$a, "is not true"))
	}
}

omxCheckCloseEnough <- function(a, b, epsilon=10^(-15)) {	
	if(abs(a - b) > epsilon) {
		stop(paste("Error:", a, "and", b, "are equal to within", epsilon))
	}
}

omxCheckWithinPercentError <- function(a, b, epsilon=10^(-15)) {	
	if(abs((a - b)/a*100) > epsilon) {
		stop(paste("Error: ", b, "does not estimate", a, "within", epsilon, "percent error"))
	}
}