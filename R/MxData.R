setClassUnion("MxData", c("data.frame", "matrix"))

omxRemoveDataAliases <- function(datalist) {
	predicate <- sapply(datalist, function(x) {
		is(x, "MxData")	
	})
	if (length(predicate) > 0) {
		return(datalist[predicate])
	} else {
		return(list())
	}
}

omxDataIndex <- function(name, datalist) {
	if (is.character(datalist[[name]])) {
		index <- datalist[[name]]	
	} else {
		index <- name
	}
	noAlias <- omxRemoveDataAliases(datalist)
	return(match(index, names(noAlias)) - 1)
}