#
#   Copyright 2007-2021 by the individuals mentioned in the source code history
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

setClass(Class = "MxFitFunctionGREML",
         slots=c(
           dV = "MxCharOrNumber",
           dVnames = "character",
           MLfit = "numeric",
           REMLfit = "numeric",
           numObsAdjust = "integer",
           aug = "MxCharOrNumber",
           augGrad = "MxCharOrNumber",
           augHess = "MxCharOrNumber",
           autoDerivType = "character",
           infoMatType = "character"),
         contains = "MxBaseFitFunction")


setMethod("initialize", "MxFitFunctionGREML",
          function(.Object, ...) {
            .Object <- callNextMethod()
            .Object@dV <- ..1
            .Object@dVnames <- as.character(names(..1))
            .Object@MLfit <- numeric(0)
            .Object@REMLfit <- numeric(0)
            .Object@vector <- FALSE
            .Object@numObsAdjust <- 0L
            .Object@aug <- ..2
            .Object@augGrad <- ..3
            .Object@augHess <- ..4
            .Object@autoDerivType <- ..5
            .Object@infoMatType <- ..6
            .Object
          })

setMethod("qualifyNames", signature("MxFitFunctionGREML"),
          function(.Object, modelname, namespace) {
            .Object <- callNextMethod()
            if(length(.Object@dV)){
              .Object@dV <- sapply(.Object@dV, imxConvertIdentifier, modelname, namespace)
              .Object@dVnames <- names(.Object@dV)
            }
            if(length(.Object@aug)){.Object@aug <- imxConvertIdentifier(.Object@aug[1],modelname,namespace)}
            if(length(.Object@augGrad)){
            	.Object@augGrad <- imxConvertIdentifier(.Object@augGrad[1],modelname,namespace)
            }
            if(length(.Object@augHess)){
            	.Object@augHess <- imxConvertIdentifier(.Object@augHess[1],modelname,namespace)
            }
            return(.Object)
          })

setMethod("genericFitRename", signature("MxFitFunctionGREML"),
          function(.Object, oldname, newname) {
            if(length(.Object@dV)){
              .Object@dV <- sapply(.Object@dV, renameReference, oldname, newname)
            }
          	if(length(.Object@aug)){.Object@aug <- renameReference(.Object@aug[1], oldname, newname)}
          	if(length(.Object@augGrad)){
          		.Object@augGrad <- renameReference(.Object@augGrad[1], oldname, newname)
          	}
          	if(length(.Object@augHess)){
          		.Object@augHess <- renameReference(.Object@augHess[1], oldname, newname)
          	}
            return(.Object)
          })

setMethod("genericFitConvertEntities", "MxFitFunctionGREML",
          function(.Object, flatModel, namespace, labelsData) {

            name <- .Object@name
            modelname <- imxReverseIdentifier(flatModel, .Object@name)[[1]]
            expectName <- paste(modelname, "expectation", sep=".")

            expectation <- flatModel@expectations[[expectName]]
            dataname <- expectation@data

            return(flatModel)
          })


setMethod("genericFitFunConvert", "MxFitFunctionGREML",
          function(.Object, flatModel, model, labelsData, dependencies) {
            .Object <- callNextMethod()
            name <- .Object@name
            modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
            expectName <- paste(modelname, "expectation", sep=".")
            if (expectName %in% names(flatModel@expectations)) {
              expectIndex <- imxLocateIndex(flatModel, expectName, name)
            } else {
              expectIndex <- as.integer(NA)
            }
            .Object@expectation <- expectIndex
            if(length(.Object@dV)){
              .Object@dV <- sapply(.Object@dV, imxLocateIndex, model=flatModel, referant=name)
            }
            if(length(.Object@aug)){.Object@aug <- imxLocateIndex(.Object@aug[1], model=flatModel, referant=name)}
            if(length(.Object@augGrad)){
            	.Object@augGrad <- imxLocateIndex(.Object@augGrad[1], model=flatModel, referant=name)
            }
            if(length(.Object@augHess)){
            	.Object@augHess <- imxLocateIndex(.Object@augHess[1], model=flatModel, referant=name)
            }
            return(.Object)
          })


setMethod("genericFitInitialMatrix", "MxFitFunctionGREML",
          function(.Object, flatModel) {return(matrix(as.double(NA), 1, 1))})

setMethod("generateReferenceModels", "MxFitFunctionGREML",
					function(.Object, model, distribution) {
						stop("Reference models for GREML expectation are not implemented")
					})


mxFitFunctionGREML <- function(
	dV=character(0), aug=character(0), augGrad=character(0), augHess=character(0), autoDerivType=c("semiAnalyt","numeric"),
	infoMatType=c("average","expected")){
	autoDerivType = as.character(match.barg(autoDerivType,c("semiAnalyt","numeric")))
	infoMatType = as.character(match.barg(infoMatType,c("average","expected")))
  return(new("MxFitFunctionGREML",dV=dV,aug=aug,augGrad=augGrad,augHess=augHess,autoDerivType=autoDerivType,infoMatType=infoMatType))
}
