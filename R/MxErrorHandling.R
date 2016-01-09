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

matchStack <- function(syscall, function_name) {
    return(identical(syscall[[1]], function_name))
}

##' imxLocateFunction
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @param function_name function_name
imxLocateFunction <- function(function_name) {
    callstack <- sys.calls()
    if (is.null(callstack)) {
        msg <- paste("(oops) could not find function",
            function_name)
        return(msg)
    }
    query <- sapply(callstack, matchStack, as.symbol(function_name))
    matches <- which(query)
    if (length(matches) == 0) {
      #If the function being sought is not actually in the call stack...
      msg <- paste("(oops) could not find function",
                   function_name)
      return(msg)
    } else {
        firstmatch <- matches[[1]]
        return(callstack[[firstmatch]])
    }
}



