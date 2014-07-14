#
#   Copyright 2007-2014 The OpenMx Project
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

mxVersion <- function (model=NULL, verbose=TRUE) {
    pvers <- try(packageVersion("OpenMx"))
    if ("try-error" %in% class(pvers)) {
        pvers = NA
    }
	if(verbose){
	    if ("try-error" %in% class(pvers)) {
			msg = paste0("OpenMx version: unknown - please report this to http://openmx.psyc.virginia.edu/forums")
	    }else{
			msg = paste0("OpenMx version: ", pvers)	    	
	    }

		msg = paste0(msg, "\nR version: ", version$version.string)	    	
		msg = paste0(msg, "\nPlatform: ", version$platform)	    	

		msg = paste0(msg, "\nDefault optimiser: ", mxOption(NULL, "Default optimizer"))
		if (!is.null(model)) {
		    msg = paste0(msg, "(optimizer for this model is ", mxOption(model, "Default optimizer"), ")")
		}
	    message(msg)
	}
	invisible(pvers)
}
