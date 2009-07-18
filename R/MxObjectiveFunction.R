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


#
# The virtual base class for all objective functions
#
setClass(Class = "MxBaseObjective", 
	representation = representation(
		name = "character",
		data = "MxCharOrNumber",
		result = "matrix", "VIRTUAL"))

setClassUnion("MxObjective", c("NULL", "MxBaseObjective"))

setGeneric("omxObjFunNamespace", 
	function(.Object, modelname, namespace) {
	return(standardGeneric("omxObjFunNamespace"))
})

setGeneric("omxObjFunConvert", 
	function(.Object, flatModel, model) {
	return(standardGeneric("omxObjFunConvert"))	
})

setGeneric("omxObjModelConvert",
	function(.Object, model) {
	return(standardGeneric("omxObjModelConvert"))
})

# NA indicates don't make any changes to the model
setMethod("omxObjModelConvert", "MxBaseObjective",
	function(.Object, model) {
		return(NA)
})