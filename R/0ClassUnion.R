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


##' A character or integer
##' @name MxCharOrNumber-class
setClassUnion("MxCharOrNumber", c("character", "integer"))

##' An optional character
##' @name MxOptionalChar-class
setClassUnion("MxOptionalChar", c("NULL", "character"))

##' @title An optional logical
##' @name MxOptionalLogical-class
##' @rdname MxOptionalLogical-class
##' @description This is an internal class, the union of NULL and logical.
setClassUnion("MxOptionalLogical", c("NULL", "logical"))

##' A character, integer, or NULL
##' @name MxOptionalCharOrNumber-class
setClassUnion("MxOptionalCharOrNumber", c("NULL", "character", "integer"))

##' An optional list
##' @name MxListOrNull-class
setClassUnion("MxListOrNull", c("list", "NULL"))

##' A character, list or NULL
##' @name MxCharOrList-class
setClassUnion("MxCharOrList", c("character", "list"))

##' An optional matrix
##' @name MxOptionalMatrix-class
setClassUnion("MxOptionalMatrix", c("NULL", "matrix"))

##' An optional numeric
##' @name MxOptionalNumeric-class
setClassUnion("MxOptionalNumeric", c("NULL", "numeric"))
