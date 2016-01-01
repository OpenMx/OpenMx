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


#
# The virtual base class for all objective functions meta-data
#
##' @title MxBaseObjectiveMetaData
##' @name MxBaseObjectiveMetaData-class
##'
##' @description
##' This is an internal class and should not be used directly.
##' It is the virtual base class for all objective functions meta-data
##'
##' @aliases
##' MxBaseObjectiveMetaData
##' @rdname MxBaseObjectiveMetaData-class
setClass(Class = "MxBaseObjectiveMetaData", 
	representation = representation("VIRTUAL"))

