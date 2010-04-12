#
#   Copyright 2007-2010 The OpenMx Project
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

mxCompare <- function(base, comparison, digits = 3) {
	if (missing(base)) {
		stop("'base' argument be a MxModel object or list of MxModel objects")	
	}
	if (missing(comparison)) {
		stop("'comparison' argument be a MxModel object or list of MxModel objects")	
	}
	if (length(digits) != 1 || is.numeric(digits) == FALSE || is.na(digits) == TRUE) {
		stop("'digits' argument must be a numeric value")
	}
	if (is.list(base)) {
		base <- unlist(base)
	} else {
		base <- list(base)
	}
	if (is.list(comparison)) {
		comparison <- unlist(comparison)
	} else {
		comparison <- list(comparison)
	}
	if(!all(sapply(base, is, "MxModel"))) {
		stop("The 'base' argument must consist of MxModel objects")
	}
	if(!all(sapply(comparison, is, "MxModel"))) {
		stop("The 'comparison' argument must consist of MxModel objects")
	}
	baseSummaries <- omxLapply(base, summary)	
	compareSummaries <- omxLapply(comparison, summary)
	resultsTable <- showFitStatistics(baseSummaries, compareSummaries, digits)
	return(resultsTable)
}

showFitStatistics <- function(baseSummaries, compareSummaries, digits) {
	statistics <- mapply(collectStatistics, baseSummaries, compareSummaries, digits)
	statistics <- data.frame(t(statistics))
	return(statistics)
}

collectStatistics <- function(refSummary, otherSummary, digits) {
	otherStats <- data.frame(stringsAsFactors = FALSE,
		refSummary$modelName,
		otherSummary$modelName,
		otherSummary$estimatedParameters,
		signif(otherSummary$Minus2LogLikelihood, digits),
		otherSummary$degreesOfFreedom,
		signif(otherSummary$AIC.Mx, digits),
		signif(otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood, digits),
		otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom, 
		signif(pchisq(otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood,
			otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom, lower.tail=FALSE), digits))
	names(otherStats) <- c("base", "comparison", "ep", "minus2LL", "df", "AIC", "diffLL", "diffdf", "p")
	return(otherStats)
}