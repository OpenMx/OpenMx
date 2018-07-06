#
#   Copyright 2007-2018 by the individuals mentioned in the source code history
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
           numObs = "integer",
           aug = "MxCharOrNumber",
           augGrad = "MxCharOrNumber",
           augHess = "MxCharOrNumber"),
         contains = "MxBaseFitFunction")


setMethod("initialize", "MxFitFunctionGREML",
          function(.Object, name = 'fitfunction', dV=character(0), MLfit=0, vector=FALSE, numObs=0L, aug=character(0),
          				 augGrad=character(0), augHess=character(0)) {
            .Object@name <- name
            .Object@dV <- dV
            .Object@dVnames <- as.character(names(dV))
            .Object@MLfit <- MLfit
            .Object@vector <- vector
            .Object@numObs <- numObs
            .Object@aug <- aug
            .Object@augGrad <- augGrad
            .Object@augHess <- augHess
            return(.Object)
          }
)


setMethod("qualifyNames", signature("MxFitFunctionGREML"), 
          function(.Object, modelname, namespace) {
            .Object@name <- imxIdentifier(modelname, .Object@name)
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


mxFitFunctionGREML <- function(dV=character(0), aug=character(0), augGrad=character(0), augHess=character(0)){
  return(new("MxFitFunctionGREML",dV=dV,aug=aug,augGrad=augGrad,augHess=augHess))
}
