#
# Objective Functions are for optimizing!
#
setClass(Class = "MxBaseObjective", 
	representation = representation(
		name = "character",
		result = "Matrix", "VIRTUAL"))

setClassUnion("MxObjective", c("NULL", "MxBaseObjective"))

setGeneric("omxObjFunConvert", function(.Object, model) {
	return(standardGeneric("omxObjFunConvert"))	
})