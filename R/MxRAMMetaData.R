#
#   Copyright 2007-2016 The OpenMx Project
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

##' @title Meta Data for RAM
##' @name MxRAMMetaData-class
##' @rdname MxRAMMetaData-class
##' @description This is an internal class, the meta data for RAM.
setClass(Class = "MxRAMMetaData",
	representation = representation(
		A = "MxCharOrNumber",
		S = "MxCharOrNumber",
		F = "MxCharOrNumber",
		M = "MxCharOrNumber",
		depth = "integer"),
	contains = "MxBaseObjectiveMetaData")

setMethod("initialize", "MxRAMMetaData",
	function(.Object, A, S, F, M, depth) {
		.Object@A <- A
		.Object@S <- S
		.Object@F <- F
		.Object@M <- M
		.Object@depth <- depth
		return(.Object)
	}
)


