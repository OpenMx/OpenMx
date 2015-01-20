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

setClass(Class = "MxExpectationGREML",
         slots = c(
           V = "MxCharOrNumber",
           X = "MxCharOrNumber",
           y = "MxCharOrNumber",
           fixedEffects = "logical",
           dims = "character",
           definitionVars = "list",
           numStats = "numeric",
           name = "character"),
         contains = "MxBaseExpectation")


setMethod("initialize", "MxExpectationGREML",
          function(.Object, V, X, y, fixedEffects,
                   data = as.integer(NA), definitionVars = list(), name = 'expectation') {
            .Object@name <- name
            .Object@V <- V
            .Object@X <- X
            .Object@y <- y
            .Object@fixedEffects <- fixedEffects
            .Object@definitionVars <- definitionVars
            .Object@data <- data
            .Object@dims <- "foo"
            return(.Object)
          }
)


setMethod("qualifyNames", signature("MxExpectationGREML"), 
          function(.Object, modelname, namespace) {
            .Object@name <- imxIdentifier(modelname, .Object@name)
            .Object@V <- imxConvertIdentifier(.Object@V, modelname, namespace)
            .Object@X <- imxConvertIdentifier(.Object@X, modelname, namespace)
            .Object@y <- sapply(.Object@y, imxConvertIdentifier, modelname, namespace)
            .Object@data <- imxConvertIdentifier(.Object@data, modelname, namespace)
            return(.Object)
          })


setMethod("genericExpDependencies", signature("MxExpectationGREML"),
          function(.Object, dependencies) {
            sources <- c(.Object@V)
            sources <- sources[!is.na(sources)]
            dependencies <- imxAddDependency(sources, .Object@name, dependencies)
            return(dependencies)
          })


setMethod("genericExpAddEntities", "MxExpectationGREML",
          function(.Object, job, flatJob, labelsData) {return(job)}
)


#setMethod("genericExpConvertEntities", "MxExpectationGREML",


setMethod("genericExpRename", signature("MxExpectationGREML"),
          function(.Object, oldname, newname) {
            .Object@X <- renameReference(.Object@X, oldname, newname)
            .Object@V <- renameReference(.Object@V, oldname, newname)
            .Object@data <- renameReference(.Object@data, oldname, newname)
            .Object@y <- sapply(.Object@y, renameReference, oldname, newname)		
            return(.Object)
          })

mxExpectationGREML <- function(V, X, y, fixedEffects=TRUE) {
  fixedEffects <- as.logical(fixedEffects)
  if (missing(V) || typeof(V) != "character") {
    stop("argument 'V' is not of type 'character' (the name of the expected covariance matrix)")
  }
  if ( missing(X) || typeof(X) != "character" )  {
    stop("argument 'X' is not of type 'character' (the name of the matrix of covariates)")
  }
  if ( missing(y) || typeof(y) != "character" ) {
    stop("argument 'y' is not of type 'character'")
  }
  return(new("MxExpectationGREML", V, X, y, fixedEffects))
}


setMethod("genericExpFunConvert", "MxExpectationGREML", 
          function(.Object, flatModel, model, labelsData, defVars, dependencies) {
            modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
            name <- .Object@name
            if(length(defVars)){stop("GREML expectation incompatible with defintion variables")}
            #There just needs to be something in the data slot, since the backend expects it:
            if(single.na(.Object@data)){.Object@data <- mxData(matrix(as.double(NA),1,1,dimnames = list(" "," ")), type="raw")}
            mxDataObject <- flatModel@datasets[[.Object@data]]
            dataName <- .Object@data
            .Object@data <- imxLocateIndex(flatModel, .Object@data, name)
            #covName <- .Object@V
            #covariance <- flatModel[[covName]]
            .Object@V <- imxLocateIndex(flatModel, .Object@V, name)
            #fixefName <- .Object@X
            .Object@X <- imxLocateIndex(flatModel, .Object@X, name)
            .Object@y <- imxLocateIndex(flatModel, .Object@y, name)
            #mxDataObject <- flatModel@datasets[[.Object@data]]
            #if (inherits(mxDataObject, "MxDataDynamic")) return(.Object)
            if (mxDataObject@type != "raw") {
              stop("GREML expectation only compatible with raw data",call.=F)
            }
            checkNumericData(mxDataObject)
            if(sum(sapply(mxDataObject@observed, is.factor))>0){
              stop("GREML expectation treats ordinal variables as continuous", call.=F)
            }
            return(.Object)
          })

