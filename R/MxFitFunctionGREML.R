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

setClass(Class = "MxFitFunctionGREML", 
         slots=c(
           dV = "MxCharOrNumber",
           dVnames = "character",
           MLfit = "numeric",
           numObs = "integer",
           Aug = "MxCharOrNumber",
           AugGrad = "MxCharOrNumber",
           AugHess = "MxCharOrNumber"),
         contains = "MxBaseFitFunction")


setMethod("initialize", "MxFitFunctionGREML",
          function(.Object, name = 'fitfunction', dV=character(0), MLfit=0, vector=FALSE, numObs=0L, Aug=character(0),
          				 AugGrad=character(0), AugHess=character(0)) {
            .Object@name <- name
            .Object@dV <- dV
            .Object@dVnames <- as.character(names(dV))
            .Object@MLfit <- MLfit
            .Object@vector <- vector
            .Object@numObs <- numObs
            .Object@Aug <- Aug
            .Object@AugGrad <- AugGrad
            .Object@AugHess <- AugHess
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
            if(length(.Object@Aug)){.Object@Aug <- imxConvertIdentifier(.Object@Aug[1],modelname,namespace)}
            if(length(.Object@AugGrad)){
            	.Object@AugGrad <- imxConvertIdentifier(.Object@AugGrad[1],modelname,namespace)
            }
            if(length(.Object@AugHess)){
            	.Object@AugHess <- imxConvertIdentifier(.Object@AugHess[1],modelname,namespace)
            }
            return(.Object)
          })

setMethod("genericFitRename", signature("MxFitFunctionGREML"),
          function(.Object, oldname, newname) {
            if(length(.Object@dV)){
              .Object@dV <- sapply(.Object@dV, renameReference, oldname, newname)
            }
          	if(length(.Object@Aug)){.Object@Aug <- renameReference(.Object@Aug[1], oldname, newname)}
          	if(length(.Object@AugGrad)){
          		.Object@AugGrad <- sapply(.Object@AugGrad[1], renameReference, oldname, newname)
          	}
          	if(length(.Object@AugHess)){
          		.Object@AugHess <- sapply(.Object@AugHess[1], renameReference, oldname, newname)
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
            if(length(.Object@Aug)){.Object@Aug <- imxLocateIndex(.Object@Aug[1], model=flatModel, referant=name)}
            if(length(.Object@AugGrad)){
            	.Object@AugGrad <- imxLocateIndex(.Object@AugGrad[1], model=flatModel, referant=name)
            }
            if(length(.Object@AugHess)){
            	.Object@AugHess <- imxLocateIndex(.Object@AugHess[1], model=flatModel, referant=name)
            }
            return(.Object)
          })


setMethod("genericFitInitialMatrix", "MxFitFunctionGREML",
          function(.Object, flatModel) {return(matrix(as.double(NA), 1, 1))})

setMethod("generateReferenceModels", "MxFitFunctionGREML",
					function(.Object, model) {
						stop("reference models for GREML expectation not implemented")
					})


mxFitFunctionGREML <- function(dV=character(0), Aug=character(0), AugGrad=character(0), AugHess=character(0)){
  return(new("MxFitFunctionGREML",dV=dV,Aug=Aug,AugGrad=AugGrad,AugHess=AugHess))
}
