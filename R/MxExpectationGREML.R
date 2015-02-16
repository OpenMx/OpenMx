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

setClass(Class = "MxExpectationGREML",
         slots = c(
           V = "MxCharOrNumber",
           X = "MxCharOrNumber",
           y = "MxCharOrNumber",
           dV = "MxCharOrNumber",
           dVnames = "character",
           numFixEff = "integer",
           dims = "character",
           definitionVars = "list",
           numStats = "numeric",
           name = "character"),
         contains = "MxBaseExpectation")

#Since there is a small chance that future developers or sophisticated users might create MxExpectationGREML
#objects with new() instead of mxExpectationGREML(), the class constructor should provide defaults for all 
#the slots...
setMethod("initialize", "MxExpectationGREML",
          function(.Object, V=character(0), X=character(0), y=character(0), dV=character(0), 
                   dVnames=character(0),
                   data = as.integer(NA), definitionVars = list(), name = 'expectation') {
            .Object@name <- name
            .Object@V <- V
            .Object@X <- X
            .Object@y <- y
            .Object@dV <- dV
            .Object@numFixEff <- integer(0)
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
            if(length(.Object@dV)){
              .Object@dV <- sapply(.Object@dV, imxConvertIdentifier, modelname, namespace)
              .Object@dVnames <- names(.Object@dV)
            }
            .Object@data <- imxConvertIdentifier(.Object@data, modelname, namespace)
            return(.Object)
          })


setMethod("genericExpDependencies", signature("MxExpectationGREML"),
          function(.Object, dependencies) {
            sources <- c(.Object@V)
            sources <- c(.Object@y)
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
            if(length(.Object@dV)){
              .Object@dV <- sapply(.Object@dV, renameReference, oldname, newname)
            }
            return(.Object)
          })

mxExpectationGREML <- function(V, X="X", y="y", dV=character(0)){
  if (missing(V) || typeof(V) != "character") {
    stop("argument 'V' is not of type 'character' (the name of the expected covariance matrix)")
  }
  if ( missing(X) || typeof(X) != "character" )  {
    stop("argument 'X' is not of type 'character' (the name of the matrix of covariates)")
  }
  if ( missing(y) || typeof(y) != "character" ) {
    stop("argument 'y' is not of type 'character'")
  }
  return(new("MxExpectationGREML", V, X, y, dV))
}


setMethod("genericExpFunConvert", "MxExpectationGREML", 
          function(.Object, flatModel, model, labelsData, defVars, dependencies) {
            modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
            name <- .Object@name
            if(length(defVars)){stop("definition variables are incompatible (and unnecessary) with GREML expectation",call.=F)}
            #There just needs to be something in the data slot, since the backend expects it:
            if(is.na(.Object@data)){
              msg <- paste("the GREML expectation function",
                           "does not have a dataset associated with it in model",
                           omxQuotes(modelname),
                           "\nthe model just needs to contain an arbitrary MxData object")
              stop(msg, call. = FALSE)
            }
            mxDataObject <- flatModel@datasets[[.Object@data]]
            checkNumericData(mxDataObject)
            if(!is.null(mxDataObject$numObs)){mxDataObject$numObs <- 0}
            #Get number of observed statistics BEFORE call to backend, so summary() can use it:
            .Object@numStats <- nrow(mxEvalByName(.Object@y,model,T))
            #Also get number of fixed effects:
            .Object@numFixEff <- ncol(mxEvalByName(.Object@X,model,T))
            #Right now, most checks on the data are unnecessary, since GREML presently ignores mxData objects...
            #if (mxDataObject@type != "raw") {
            #  stop("GREML expectation only compatible with raw data",call.=F)
            #}
            #if(sum(sapply(mxDataObject@observed, is.factor))>0){
            #  stop("GREML expectation not compatible with ordinal data", call.=F)
            #}
            dataName <- .Object@data
            .Object@data <- imxLocateIndex(flatModel, .Object@data, name)
            .Object@V <- imxLocateIndex(flatModel, .Object@V, name)
            .Object@X <- imxLocateIndex(flatModel, .Object@X, name)
            .Object@y <- imxLocateIndex(flatModel, .Object@y, name)
            if(length(.Object@dV)){
              .Object@dV <- sapply(.Object@dV, imxLocateIndex, model=flatModel, referant=name)
            }
            return(.Object)
          })

