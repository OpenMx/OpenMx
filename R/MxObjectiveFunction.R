#
# Objective Functions are for optimizing!
#
setClass(Class = "MxObjective", 
	representation = representation(
		name = "character", "VIRTUAL"))

setClassUnion("MxCharOrNumber", c("character", "numeric"))

setGeneric("omxObjFunConvert", function(.Object, model) {
	return(standardGeneric("omxObjFunConvert"))	
})