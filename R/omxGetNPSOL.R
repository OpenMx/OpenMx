#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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

# nocov start

##' omxGetNPSOL
##'
##' Get the non-CRAN version of OpenMx from the OpenMx website.
##' 
##' @details
##' This function 
##' @return
##' Invisible NULL
##' @examples
##' \donttest{omxGetNPSOL()}
##'

omxGetNPSOL <- function() {
	if (imxHasNPSOL()) {
		message(paste("NPSOL is available.",
			      "You have already installed the non-CRAN version of OpenMx."))
		return()
	}
	if(.Platform$OS.type=="windows"){
		message(
			paste("Windows users must clear R's workspace, restart R, run\n",
						"source('http://openmx.ssri.psu.edu/getOpenMx.R')\n",
						", and then restart R again before trying to load OpenMx.")
		)
		return()
	}
    if(version$major < 3) {
        message(paste0("You are using R 2.15 or earlier.  ",
            "OpenMx 2.0 and higher do not support versions ",
            "of R before 3.0, so I'm fetching OpenMx 1.4 instead.\n",
            "Getting OpenMx 1.4 from http://openmx.ssri.psu.edu/."))
        source("http://openmx.ssri.psu.edu/getOpenMx1.4.R")
    } else {
        message("Getting OpenMx from http://openmx.ssri.psu.edu/.")
        source("http://openmx.ssri.psu.edu/getOpenMx.R")
    }
}
# nocov end
