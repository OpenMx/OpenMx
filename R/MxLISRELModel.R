#
#   Copyright 2007-2013 The OpenMx Project
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

setClass(Class = "MxLISRELModel",
	representation = representation(),
	contains = "MxModel")

imxModelTypes[['LISREL']] <- "MxLISRELModel"

# imxVariableTypes <- c(imxVariableTypes, "exogenous", "endogenous")


#--------------------------------------------------------------------
# Define generic functions

setMethod("imxTypeName", "MxLISRELModel", 
	function(model) { "LISREL" }
)

setMethod("imxInitModel", "MxLISRELModel", 
	function(model) {
		stop("Not implemented")
		
		if (is.null(model[['LX']])) {
			model[['LX']] <- createMatrixLoa(model, TRUE)
		}
		if (is.null(model[['LY']])) {
			model[['LY']] <- createMatrixLoa(model, FALSE)
		}
		if (is.null(model[['BE']])) {
			model[['BE']] <- createMatrixReg(model, FALSE)
		}
		if (is.null(model[['GA']])) {
			model[['GA']] <- createMatrixReg(model, TRUE)
		}
		if (is.null(model[['PH']])) {
			model[['PH']] <- createMatrixErr(model, TRUE)
		}
		if (is.null(model[['PS']])) {
			model[['PS']] <- createMatrixErr(model, FALSE)
		}
		if (is.null(model[['TD']])) {
			model[['TD']] <- createMatrixRes(model, TRUE)
		}
		if (is.null(model[['TE']])) {
			model[['TE']] <- createMatrixRes(model, FALSE)
		}
		if (is.null(model[['TH']])) {
			model[['TH']] <- createMatrixRes(model, 'both')
		}
		if (is.null(model[['TX']])) {
			model[['TX']] <- createMatrixInt(model, TRUE)
		}
		if (is.null(model[['TY']])) {
			model[['TY']] <- createMatrixInt(model, FALSE)
		}
		if (is.null(model[['KA']])) {
			model[['KA']] <- createMatrixMea(model, TRUE)
		}
		if (is.null(model[['AL']])) {
			model[['AL']] <- createMatrixMea(model, FALSE)
		}
		if (is.null(model[['expectation']])) {
			model[['expectation']] <- mxExpectationLISREL('LX', 'LY', 'BE', 'GA', 'PH', 'PS', 'TD', 'TE', 'TH', 'TX', 'TY', 'KA', 'AL')
		}
		if (is.null(model[['fitfunction']])) {
			model[['fitfunction']] <- mxFitFunctionML()
		}
		return(model)
	}
)

createMatrixLoa <- function(model, exogenous) {}
createMatrixReg <- function(model, exogenous) {}
createMatrixErr <- function(model, exogenous) {}
createMatrixRes <- function(model, exogenous) {}
createMatrixMea <- function(model, exogenous) {}
createMatrixInt <- function(model, exogenous) {}

setMethod("imxModelBuilder", "MxLISRELModel", 
	function(model, lst, name, 
		manifestVars, latentVars, lst.call, remove, independent) {
		stop("Not implemented")
	}
)

setMethod("imxVerifyModel", "MxLISRELModel",
	function(model) {
		return(TRUE)
	}
)


setReplaceMethod("[[", "MxLISRELModel",
	function(x, i, j, value) {
		stop("Not implemented")
	}
)

setReplaceMethod("$", "MxLISRELModel",
	function(x, name, value) {
		stop("Not implemented")
	}
)
