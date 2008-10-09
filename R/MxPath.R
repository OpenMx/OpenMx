mxCreatePath <- function(from, to, all = FALSE, free = NULL, 
	arrows = NULL, startVal = NULL, 
	endVal = NULL, algebra = NULL,
	name = NULL, label = NULL,
	boundMax = NULL, boundMin = NULL,
	ciUpper = NULL, ciLower = NULL) {
	if(length(from) == 1 && length(to) == 1) {
		return(mxSinglePath(from, to, free, 
			arrows, startVal, endVal, 
			algebra, name, label, boundMax, 
			boundMin, ciUpper, ciLower))
	} else {
		from <- as.list(from)
		to <- as.list(to)
		if (all) {
			from <- rep(from, each=length(to))	
		}
		result <- suppressWarnings(mapply(mxSinglePath, from, to,
			free, arrows, startVal, endVal,
				algebra, name, label, boundMax,
				boundMin, ciUpper, ciLower, SIMPLIFY=FALSE))
		return(result)			
	}
}

mxSinglePath <- function(from, to, free = NULL, 
	arrows = NULL, startVal = NULL, 
	endVal = NULL, algebra = NULL,
	name = NULL, label = NULL,
	boundMax = NULL, boundMin = NULL,
	ciUpper = NULL, ciLower = NULL) {
	result <- list()
	result[['from']] <- from
	result[['to']] <- to
	result[['free']] <- free
	result[['arrows']] <- arrows
	result[['startVal']] <- startVal
	result[['endVal']] <- endVal
	result[['algebra']] <- algebra
	result[['name']] <- name
	result[['label']] <- label		
	result[['boundMax']] <- boundMax
	result[['boundMin']] <- boundMin	
	result[['ciUpper']] <- ciUpper
	result[['ciLower']] <- ciLower
	return(result)
}

isMxPath <- function(value) {
	return(is.list(value) && 
		!is.null(value[['from']]) &&
		!is.null(value[['to']]))
}