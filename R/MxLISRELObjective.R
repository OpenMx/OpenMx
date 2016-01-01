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

#--------------------------------------------------------------------
# Author: Michael D. Hunter
# Filename: MxLISRELObjective.R
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Revision History
#   Mon Feb 20 13:03:21 Central Standard Time 2012 -- Michael Hunter added means
#   Sat Apr 07 19:48:33 Central Daylight Time 2012 -- Michael Hunter added lots of error checking
# 

#--------------------------------------------------------------------
# **DONE**
mxLISRELObjective <- function(LX=NA, LY=NA, BE=NA, GA=NA, PH=NA, PS=NA, TD=NA, TE=NA, TH=NA, TX = NA, TY = NA, KA = NA, AL = NA, dimnames = NA, thresholds = NA, vector = FALSE, threshnames = dimnames) {
	LX <- checkLISRELargument(LX, "LX")
	LY <- checkLISRELargument(LY, "LY")
	BE <- checkLISRELargument(BE, "BE")
	GA <- checkLISRELargument(GA, "GA")
	PH <- checkLISRELargument(PH, "PH")
	PS <- checkLISRELargument(PS, "PS")
	TD <- checkLISRELargument(TD, "TD")
	TE <- checkLISRELargument(TE, "TE")
	TH <- checkLISRELargument(TH, "TH")
	TX <- checkLISRELargument(TX, "TX")
	TY <- checkLISRELargument(TY, "TY")
	KA <- checkLISRELargument(KA, "KA")
	AL <- checkLISRELargument(AL, "AL")
	
	if (single.na(thresholds)) thresholds <- as.character(NA)
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("Dimnames argument is not a character vector")
	}
	if (length(thresholds) != 1) {
		stop("Thresholds argument must be a single matrix or algebra name")
	}
	if (length(dimnames) == 0) {
		stop("Dimnames argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for dimnames vector")
	}
	if (length(vector) > 1 || typeof(vector) != "logical") {
		stop("Vector argument is not a logical value")
	}
	threshnames <- checkThreshnames(threshnames)
	expectation <- mxExpectationLISREL(LX, LY, BE, GA, PH, PS, TD, TE, TH, TX, TY, KA, AL, dimnames, thresholds, threshnames)
	fitfunction <- mxFitFunctionML(vector)
	msg <- paste("Objective functions have been deprecated.",
		"Please use mxExpectationLISREL() and mxFitFunctionML() instead.")
	warning(msg)
	return(list(expectation=expectation, fitfunction=fitfunction))
}



