#
#   Copyright 2007-2012 The OpenMx Project
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
# 


#-------------------------------------------------------------------------------------
# TODO Add a data argument to allow users to give mxSaturatedModel a data set
# TODO Add ability to fit models to ordinal data
# TODO Add ability to fit multiple groups


#-------------------------------------------------------------------------------------
# Saturated Model function definition

omxSaturatedModel <- function(model) {
	if (!(isS4(model) && is(model, "MxModel"))) {
		stop("'model' argument must be a MxModel object")
	}
	datasource <- model$data
	if (is.null(datasource)) {
		stop("'model' argument does not contain any data")
	}
	numVar <- ncol(datasource@observed)
	varnam <- colnames(datasource@observed)
	if(is.null(varnam)) {
		varnam <- paste("x", 1:numVar, sep="")
		dimnames(datasource@observed) <- list(varnam, varnam)
	}
	if(datasource@type == "raw") {
		startcov <- t(chol(cov(datasource@observed, use="pairwise.complete.obs")))
		startcov <- startcov[lower.tri(startcov, TRUE)]
		startmea <- colMeans(datasource@observed, na.rm=TRUE)
	} else {
		startcov <- 0.3
		startmea <- 3.0
	}
	saturatedModel <- mxModel(
		name=paste("Saturated", model@name),
		datasource,
		mxMatrix(type="Lower",
			nrow=numVar,
			ncol=numVar,
			values=startcov,
			free=T,
			name="ltCov"),
		mxAlgebra(name="satCov", expression= ltCov %*% t(ltCov), dimnames=list(varnam, varnam))
	)
	if(datasource@type == "raw" || !any(is.na(datasource@means)) ) {
		saturatedModel <- mxModel(saturatedModel,
			mxMatrix(nrow=1, ncol=numVar, values=startmea, free=T, name="satMea", dimnames=list(NA, varnam)),
			mxMLObjective("satCov", "satMea")
		)
	} else {
		saturatedModel <- mxModel(saturatedModel,
			mxMLObjective("satCov")
		)
	}
	
	saturatedModel <- mxOption(saturatedModel, "Calculate Hessian", "No")
	saturatedModel <- mxOption(saturatedModel, "Standard Errors", "No")
	return(saturatedModel)
}


#-------------------------------------------------------------------------------------
# End

