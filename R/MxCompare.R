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

mxCompare <- function(base, comparison, ..., transform = signif, digits = 3, all = FALSE) {
	if (missing(base)) {
		stop("'base' argument be a MxModel object or list of MxModel objects")	
	}
	if (missing(comparison)) {
		stop("'comparison' argument be a MxModel object or list of MxModel objects")	
	}
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxCompare does not accept values for the '...' argument")
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
	resultsTable <- showFitStatistics(baseSummaries, compareSummaries, transform, digits, all)
	return(resultsTable)
}

showFitStatistics <- function(baseSummaries, compareSummaries, transform, digits, all) {
	statistics <- list()
	if(all) {
		for(i in 1:length(baseSummaries)) {
			nextBaseSummary <- baseSummaries[[i]]
			statistics[[length(statistics) + 1]] <- collectBaseStatistics(nextBaseSummary, 
				transform, digits)
			for(j in 1:length(compareSummaries)) {
				nextCompareSummary <- compareSummaries[[j]]
				statistics[[length(statistics) + 1]] <- collectStatistics(nextBaseSummary,
					nextCompareSummary, transform, digits)
			}
		}
	} else {
		maxLength <- max(length(baseSummaries), length(compareSummaries))
		previousBaseSummaryIndex <- -1
		for(i in 1:maxLength) {
			nextBaseSummaryIndex <- (i - 1) %% length(baseSummaries) + 1
			nextCompareSummaryIndex <- (i - 1) %% length(compareSummaries) + 1
			nextBaseSummary <- baseSummaries[[nextBaseSummaryIndex]]
			nextCompareSummary <- compareSummaries[[nextCompareSummaryIndex]]
			if (previousBaseSummaryIndex != nextBaseSummaryIndex) {
				statistics[[length(statistics) + 1]] <- collectBaseStatistics(nextBaseSummary,
					transform, digits)
			}
			statistics[[length(statistics) + 1]] <- collectStatistics(nextBaseSummary,
				nextCompareSummary, transform, digits)
			previousBaseSummaryIndex <- nextBaseSummaryIndex
		}
	}
	statistics <- do.call(rbind, statistics)
	return(statistics)
}

collectBaseStatistics <- function(refSummary, transform, digits) {
	baseStats <- data.frame(stringsAsFactors = FALSE,
		refSummary$modelName,
		as.character(NA),
		refSummary$estimatedParameters,
		transform(refSummary$Minus2LogLikelihood, digits),
		refSummary$degreesOfFreedom,
		transform(refSummary$AIC.Mx, digits),
		as.numeric(NA),
		as.numeric(NA), 
		as.numeric(NA))
	names(baseStats) <- c("base", "comparison", "ep", "minus2LL", "df", "AIC", "diffLL", "diffdf", "p")
	return(baseStats)
}

collectStatistics <- function(refSummary, otherSummary, transform, digits) {
	otherStats <- data.frame(stringsAsFactors = FALSE,
		refSummary$modelName,
		otherSummary$modelName,
		otherSummary$estimatedParameters,
		transform(otherSummary$Minus2LogLikelihood, digits),
		otherSummary$degreesOfFreedom,
		transform(otherSummary$AIC.Mx, digits),
		transform(otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood, digits),
		otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom, 
		transform(pchisq(otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood,
			otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom, lower.tail=FALSE), digits))
	names(otherStats) <- c("base", "comparison", "ep", "minus2LL", "df", "AIC", "diffLL", "diffdf", "p")
	return(otherStats)
}
