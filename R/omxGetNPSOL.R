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

##' omxGetNPSOL
##'
##' Get the non-CRAN version of OpenMx from the OpenMx website.
##' 
##' @details
##' This function 
##' @return
##' Invisible NULL
##' @examples
##' \dontrun{omxGetNPSOL()}
##'

omxGetNPSOL <- function() {
	if (imxHasNPSOL()) {
		message(paste("NPSOL is available.",
			      "You have already installed the non-CRAN version of OpenMx."))
		return()
	}
	if(.Platform$OS.type=="windows"){
		message(
			paste("Windows users should either restart R or run\n",
						"detach('package:OpenMx',unload=TRUE)\n",
						", and then run\n",
						"source('http://openmx.psyc.virginia.edu/getOpenMx.R')\n")
		)
		return()
	}
    if(version$major < 3) {
        message(paste0("You are using R 2.15 or earlier.  ",
            "OpenMx 2.0 and higher do not support versions ",
            "of R before 3.0, so I'm fetching OpenMx 1.4 instead.\n",
            "Getting OpenMx 1.4 from http://openmx.psyc.virginia.edu/."))
        source("http://openmx.psyc.virginia.edu/getOpenMx1.4.R")
    } else {
        message("Getting OpenMx from http://openmx.psyc.virginia.edu/.")
        source("http://openmx.psyc.virginia.edu/getOpenMx.R")
    }
}
