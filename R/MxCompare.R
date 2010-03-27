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

mxCompare <- function(model, ..., digits = 3) {
	if(missing(model) || !is(model, "MxModel")) {
		stop("'model' argument be be a MxModel object")
	}
	otherModels <- list(...)
	otherModels <- unlist(otherModels)
	if(!all(sapply(otherModels, is, "MxModel"))) {
		stop("The '...' arguments must all be MxModel objects")
	}
	resultsTable <- showFitStatistics(model, otherModels, digits)
	rownames(resultsTable) <- sapply(1:nrow(resultsTable), compareGenerateName)
	return(resultsTable)
}

showFitStatistics <- function(reference, otherModels, digits) {
	refSummary <- summary(reference)
	base <- collectBaseStatistics(refSummary, reference@name, digits)
	if (length(otherModels) == 0) {
		return(base)
	} else {
		base <- cbind(base, "diffLL" = NA, "diffdf" = NA, "p" = NA)
		others <- lapply(otherModels, collectOtherStatistics, refSummary, digits)
		others <- do.call(rbind, others)
		retval <- rbind(base, others)
		return(retval)
	}
}

collectBaseStatistics <- function(refSummary, refName, digits) {
	stats <- data.frame(
			refName,
			refSummary$estimatedParameters,
			signif(refSummary$Minus2LogLikelihood, digits),
			refSummary$degreesOfFreedom,
			signif(refSummary$AIC.Mx, digits))
	names(stats) <- c("Name","ep","-2LL", "df", "AIC")
	return(stats)
}

collectOtherStatistics <- function(otherModel, refSummary, digits) {
	otherSummary <- summary(otherModel)
	otherStats <- data.frame(
		otherModel@name, 
		otherSummary$estimatedParameters,
		signif(otherSummary$Minus2LogLikelihood, digits),
		otherSummary$degreesOfFreedom,
		signif(otherSummary$AIC.Mx, digits),
		signif(otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood, digits),
		otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom, 
		signif(pchisq(otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood,
			otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom, lower.tail=FALSE), digits))
	names(otherStats) <- c("Name","ep","-2LL", "df", "AIC","diffLL","diffdf","p")
	return(otherStats)
}


compareGenerateName <- function(index) {
	paste("Model", index, ":")
}
