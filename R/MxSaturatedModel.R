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


#-------------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2012.10.23
# Filename: mxSaturatedModel.R
# Purpose: This is a helper function for making and fitting a saturated model.
#  It takes an mxModel and gives back the fitted saturated model.
#  Example usage
#    amod <- mxModel(blah)
#    smod <- mxSaturatedModel(amod)
#    summary(amod, SaturatedLikelihood=smod)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# Revision History
# Tue Oct 23 00:55:48 Central Daylight Time 2012 -- Michael Hunter copied file from personal work
# Wed 19 Jun 2013 12:23:45 Central Daylight Time -- Michael Hunter added ordinal support etc with aid of Mike Neale and Ryne Estabrook
# 


#-------------------------------------------------------------------------------------
# 1. Add ability to fit multiple groups
# 2. If given a model whose data set has 10 variables but the model only uses 2
#  variables, adjust function to make a saturated model of only the 2 used variables.
# TODO Check that the above are done reasonably correct.
# TODO Improve interaction between 1 & 2.
# Added ability to do Independence model.



#-------------------------------------------------------------------------------------
# Saturated Model function definition

generateNormalReferenceModels <- function(modelName, obsdata, datatype, withMeans=FALSE) {
	datasource <- mxData(observed=obsdata, type=datatype)
	numVar <- ncol(obsdata)
	varnam <- colnames(obsdata)
	if(is.null(varnam)) {
		varnam <- paste("V", 1:numVar, sep="")
		dimnames(obsdata) <- list(varnam, varnam)
	}
	if(datatype == "raw") {
		if (is.data.frame(obsdata)) {
			ordinalCols <- sapply(obsdata, is.ordered)
		}
		if(!any(ordinalCols)){
			sampcov <- cov(obsdata, use="pairwise.complete.obs")
			startcov <- t(chol(sampcov))
			startcov <- startcov[lower.tri(startcov, TRUE)]
			indepcov <- diag(sampcov)
			startmea <- colMeans(obsdata, na.rm=TRUE)
		}
		else {
			ordinalLevels <- lapply(obsdata[,ordinalCols], levels)
			numOrdinal <- sum(ordinalCols)
			maxLevels <- max(sapply(ordinalLevels, length))
			numThresholds <- maxLevels-1
			startcov <- t(chol(diag(1, numVar)))
			startcov <- startcov[lower.tri(startcov, TRUE)]
			indepcov <- diag(1, numVar)
			startmea <- rep(0, numVar)
		}
	} else {
		startcov <- 0.3
		indepcov <- 0.3
		startmea <- 3.0
	}

	ltCov <- mxMatrix(type="Lower", nrow=numVar, ncol=numVar,
			  values=startcov, free=TRUE, name="ltCov")
	saturatedModel <- mxModel(name=paste("Saturated", modelName),
				  datasource,
				  ltCov,
				  mxAlgebra(name="satCov", expression= ltCov %*% t(ltCov), dimnames=list(varnam, varnam)),
				  mxExpectationNormal("satCov"),
				  mxFitFunctionML())

	indCov <- mxMatrix(type="Diag", nrow=numVar, ncol=numVar, values=indepcov, free=TRUE,
			   name="indCov", dimnames=list(varnam, varnam))
	independenceModel <- mxModel(name=paste("Independence", modelName),
				     datasource, indCov,
				     mxExpectationNormal("indCov"), mxFitFunctionML())

	if(datatype == "raw" || withMeans) {
		saturatedModel <- mxModel(saturatedModel,
			mxMatrix(nrow=1, ncol=numVar, values=startmea, free=TRUE, name="satMea", dimnames=list(NA, varnam)),
			mxExpectationNormal("satCov", "satMea")
		)
		independenceModel <- mxModel(independenceModel,
			mxMatrix(nrow=1, ncol=numVar, values=startmea, free=TRUE, name="satMea", dimnames=list(NA, varnam)),
			mxExpectationNormal("indCov", "satMea")
		)
		if(any(ordinalCols)) {
			unitLower <- mxMatrix("Lower", numThresholds, numThresholds, values=1, free=FALSE, name="unitLower")
			thresholdDeviations <- mxMatrix("Full", 
					name="thresholdDeviations", nrow=numThresholds, ncol=numOrdinal,
					values=.2,
					free = TRUE, 
					lbound = rep( c(-Inf,rep(.01, (numThresholds-1))) , numOrdinal), # TODO adjust increment value
					dimnames = list(c(), varnam[ordinalCols]), # TODO Add threshold names
							)
			saturatedModel <- mxModel(saturatedModel,
				mxMatrix(nrow=1, ncol=numVar, values=startmea, free=c(!ordinalCols), name="satMea", dimnames=list(NA, varnam)),
				thresholdDeviations, unitLower,
				mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
				mxExpectationNormal("satCov", "satMea", thresholds="thresholdMatrix")
			)
			independenceModel <- mxModel(independenceModel,
				mxMatrix(nrow=1, ncol=numVar, values=startmea, free=c(!ordinalCols), name="satMea", dimnames=list(NA, varnam)),
				thresholdDeviations, unitLower,
				mxAlgebra(unitLower %*% thresholdDeviations, name="thresholdMatrix"),
				mxExpectationNormal("indCov", "satMea", thresholds="thresholdMatrix")
			)
		}
	}
	return(list(Saturated=saturatedModel, Independence=independenceModel))
}

ReferenceModelHelper <- function(x) {
	if ( (!(isS4(x) && is(x, "MxModel"))) && !is.data.frame(x) && !(is.matrix(x) && is.numeric(x)) ) {
		stop("The 'x' argument must be (1) an MxModel object, (2) a raw data frame, or (3) a raw data matrix.")
	}
	if ( is(x, "MxModel") ) {
		if (is.null(x$fitfunction)) {
			stop("Model", omxQuotes(x$name), "has no fitfunction")
		}
		generateReferenceModels(x$fitfunction, x)
	} else {
		obsdata <- x
		if(ncol(obsdata) != nrow(obsdata)) {
			datatype <- "raw"
		}
		else {datatype <- "cov"}
		generateNormalReferenceModels("Data Model", obsdata, datatype)
	}
}

omxSaturatedModel <- function(x, run=FALSE) {
	models <- lapply(ReferenceModelHelper(x), function(model) {
		if (!isS4(model)) return(model)
		model <- mxOption(model, "Standard Errors", "No")
		model <- mxOption(model, "Calculate Hessian", "No")
		if (run) {
			model <- mxRun(model)
		}
		model
	})
	models
}
