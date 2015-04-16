#
#   Copyright 2007-2015 The OpenMx Project
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
           numObs = "integer"),
         contains = "MxBaseFitFunction")


setMethod("initialize", "MxFitFunctionGREML",
          function(.Object, name = 'fitfunction', dV=character(0)) {
            .Object@name <- name
            .Object@dV <- dV
            .Object@vector <- FALSE
            .Object@numObs <- 0L
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
            return(.Object)
          })

setMethod("genericFitRename", signature("MxFitFunctionGREML"),
          function(.Object, oldname, newname) {
            if(length(.Object@dV)){
              .Object@dV <- sapply(.Object@dV, renameReference, oldname, newname)
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
          function(.Object, flatModel, model, labelsData, defVars, dependencies) {
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
            return(.Object)
          })


setMethod("genericFitInitialMatrix", "MxFitFunctionGREML",
          function(.Object, flatModel) {return(matrix(as.double(NA), 1, 1))})


mxFitFunctionGREML <- function(dV=character(0)){
  return(new("MxFitFunctionGREML",dV=dV))
}
