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

setClassUnion("MxMaybeFunction", c("NULL", "function"))

setClass(Class = "MxReservedName",
	representation = representation(
		name = "character",
		search = "MxMaybeFunction",
		replace = "MxMaybeFunction"))

setMethod("initialize", "MxReservedName",
	function(.Object, name = character(),
		search = function(model) { return(NULL) },
		replace = function(model, value) { return(model) }) {
		.Object@name <- name
		.Object@search <- search
		.Object@replace <- replace
		return(.Object)
	}
)

omxReservedNames <- list()

omxReservedNames[['data']] <- new("MxReservedName", "data", NULL, NULL)

omxReservedNames[['objective']] <- new("MxReservedName", "objective", NULL, NULL)

omxReservedNames[['likelihood']] <- new("MxReservedName", "likelihood")
