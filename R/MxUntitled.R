#
#   Copyright 2007-2016 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


##' imxUntitledNumberReset
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @details Resets the \link{imxUntitledNumber} counter
imxUntitledNumberReset <- function() .Call(untitledNumberReset)

##' imxUntitledNumber
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @details Increments the untitled number counter and returns its value
imxUntitledNumber <- function() .Call(untitledNumber)

##' imxUntitledName
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @details Returns a character, the name of the next untitled entity
imxUntitledName <- function() {
	name <- paste("untitled", format(imxUntitledNumber(), scientific = FALSE),
		sep="")
	return(name)
}
