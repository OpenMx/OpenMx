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

mxCompare <- function(base, comparison, ..., all = FALSE) {
	if (missing(base)) {
		stop("'base' argument be a MxModel object or list of MxModel objects")	
	}
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxCompare does not accept values for the '...' argument")
	}
	if (is.list(base)) {
		base <- unlist(base)
	} else {
		base <- list(base)
	}
	if(!all(sapply(base, is, "MxModel"))) {
		stop("The 'base' argument must consist of MxModel objects")
	}
	baseSummaries <- omxLapply(base, summary)
	
	if(missing(comparison)) {
		# no comparison models, just make a dummy list to feed to showFitStatistics
		compareSummaries <- list()
	}else{
		if (is.list(comparison)) {
			comparison <- unlist(comparison)
		} else {
			comparison <- list(comparison)
		}
		if(!all(sapply(comparison, is, "MxModel"))) {
			stop("The 'comparison' argument must consist of MxModel objects")
		} else {
			compareSummaries <- omxLapply(comparison, summary)
		}
	}
	resultsTable <- showFitStatistics(baseSummaries, compareSummaries, all)
	return(resultsTable)
}


showFitStatistics <- function(baseSummaries, compareSummaries, all)  {
    statistics <- list()
    if (all) {
        for (i in seq_along(baseSummaries)) {
            nextBaseSummary <- baseSummaries[[i]]
            statistics[[length(statistics) + 1]] <- collectBaseStatistics(nextBaseSummary)
            for (j in seq_along(compareSummaries)) {
                nextCompareSummary <- compareSummaries[[j]]
                statistics[[length(statistics) + 1]] <- collectStatistics(nextBaseSummary, nextCompareSummary)
            }
        }
    }
    else {
		if(length(compareSummaries)==0){
	        for (i in seq_along(baseSummaries)) {
	            statistics[[length(statistics) + 1]] <- collectBaseStatistics(baseSummaries[[i]])
	        }
		} else {
	        maxLength <- max(length(baseSummaries), length(compareSummaries))
	        previousBaseSummaryIndex <- -1
	        for (i in 1:maxLength) {
	            nextBaseSummaryIndex <- (i - 1)%%length(baseSummaries) + 1
	            nextCompareSummaryIndex <- (i - 1)%%length(compareSummaries) + 1
	            nextBaseSummary <- baseSummaries[[nextBaseSummaryIndex]]
	            nextCompareSummary <- compareSummaries[[nextCompareSummaryIndex]]
	            if (previousBaseSummaryIndex != nextBaseSummaryIndex) {
	                statistics[[length(statistics) + 1]] <- collectBaseStatistics(nextBaseSummary)
	            }
	            statistics[[length(statistics) + 1]] <- collectStatistics(nextBaseSummary, nextCompareSummary)
	            previousBaseSummaryIndex <- nextBaseSummaryIndex
	        }
		}		
    }
    statistics <- do.call(rbind, statistics)
    return(statistics)
}


collectBaseStatistics <- function(refSummary) {
	baseStats <- data.frame(stringsAsFactors = FALSE,
		refSummary$modelName,
		as.character(NA),
		refSummary$estimatedParameters,
		refSummary$Minus2LogLikelihood,
		refSummary$degreesOfFreedom,
		refSummary$AIC.Mx,
		as.numeric(NA),
		as.numeric(NA), 
		as.numeric(NA))
	names(baseStats) <- c("base", "comparison", "ep", "minus2LL", "df", "AIC", "diffLL", "diffdf", "p")
	return(baseStats)
}

collectStatistics <- function(refSummary, otherSummary) {
	otherStats <- data.frame(stringsAsFactors = FALSE,
		refSummary$modelName,
		otherSummary$modelName,
		otherSummary$estimatedParameters,
		otherSummary$Minus2LogLikelihood,
		otherSummary$degreesOfFreedom,
		otherSummary$AIC.Mx,
		otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood,
		otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom, 
		pchisq(otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood,
			otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom, lower.tail=FALSE))
	names(otherStats) <- c("base", "comparison", "ep", "minus2LL", "df", "AIC", "diffLL", "diffdf", "p")
	return(otherStats)
}
