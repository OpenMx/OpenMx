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

##' omxRAMtoML
##'
##' Convert a RAM model to an ML model
##' 
##' @param model the MxModel
##' @details
##' This is a legacy function that was once used to convert RAM models to ML models
##' in the old (1.0 release of OpenMx) objective function style.
##' @return
##' an ML model with an ML objective
omxRAMtoML <- function(model) {
	namespace <- imxGenerateNamespace(model)
	job <- shareData(model)
	return(RAMtoMLHelper(model, job, namespace))
}

RAMtoMLHelper <- function(model, job, namespace) {
	model <- convertRAMtoMLModel(model, job, namespace)
	if (length(model@submodels) > 0) {
		model@submodels <- lapply(model@submodels, RAMtoMLHelper, job, namespace)
	}
	return(model)
}

createNewName <- function(model, namespace, suggestedName) {
    if (availableName(model, namespace, suggestedName)) {
		return(suggestedName)
    } else {
		return(paste(suggestedName, imxUntitledName(), sep = '_'))
    }
}

convertRAMtoMLModel <- function(model, job, namespace) {
	objective <- model$objective
	if (is.null(objective)) {
		return(model)
	}
	if (!is(objective, "MxRAMObjective")) {
		return(model)
	}
	modelname <- model@name
	aName <- imxConvertIdentifier(objective@A, modelname, namespace)
	sName <- imxConvertIdentifier(objective@S, modelname, namespace)
	fName <- imxConvertIdentifier(objective@F, modelname, namespace)
	mName <- imxConvertIdentifier(objective@M, modelname, namespace)	
	iName <- createNewName(model, namespace, 'I')
    model <- mxModel(model, mxMatrix(type="Iden", 
		nrow = nrow(job[[aName]]), name = iName))
	zName <- createNewName(model, namespace, 'Z')
	zFormula <- substitute(solve(I - A),
		list(I = as.symbol(iName), A = as.symbol(aName)))
	algebra <- eval(substitute(mxAlgebra(x, y), list(x = zFormula, y = zName)))
    model <- mxModel(model, algebra)
	covName <- createNewName(model, namespace, 'covariance')
    covFormula <- substitute(F %*% Z %*% S %*% t(Z) %*% t(F),
        list(F = as.symbol(fName), Z = as.symbol(zName),
            S = as.symbol(sName)))
    algebra <- eval(substitute(mxAlgebra(x, y),
        list(x = covFormula, y = covName)))
    model <- mxModel(model, algebra)
	
	if(!single.na(mName)) {
		meansFormula <- substitute(t(F %*% Z %*% t(M)),
			list(F = as.symbol(fName), Z = as.symbol(zName),
				M = as.symbol(mName)))
		meansName <- createNewName(model, namespace, 'means')
		algebra <- eval(substitute(mxAlgebra(x, y),
			list(x = meansFormula, y = meansName)))
		model <- mxModel(model, algebra)
	}

    translatedNames <- fMatrixTranslateNames(job[[fName]]@values, modelname)
	dataset <- job[[modelname]]@data
	if (dataset@type == 'raw') {
		objectiveType <- as.symbol('mxFIMLObjective')
		if (is.na(mName)) {
			objective <- eval(substitute(obj(covariance = x, 
				thresholds = z, vector = w, dimnames = u),
				list(x = covName, z = objective@thresholds, 
					w = objective@vector, u = translatedNames, obj = objectiveType)))
		} else {
			objective <- eval(substitute(obj(covariance = x, 
				means = y, thresholds = z, vector = w, dimnames = u),
				list(x = covName, y = meansName, z = objective@thresholds, 
					w = objective@vector, u = translatedNames, obj = objectiveType)))
		}
	} else {
		objectiveType <- as.symbol('mxMLObjective')		
		if (is.na(mName)) {
			objective <- eval(substitute(obj(covariance = x, 
				thresholds = z, dimnames = u),
				list(x = covName, z = objective@thresholds, 
					u = translatedNames, obj = objectiveType)))
		} else {
			objective <- eval(substitute(obj(covariance = x, 
				means = y, thresholds = z, dimnames = u),
				list(x = covName, y = meansName, z = objective@thresholds, 
					u = translatedNames, obj = objectiveType)))
		}
	}
    model@objective <- objective
    class(model) <- 'MxModel'
	return(model)
}
