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

##' omxConstrainMLThresholds
##'
##' Add constraint to ML model to keep thresholds in order
##' 
##' @param model the MxModel to which constraints should be added
##' @param dist unused
##' @details
##' This function adds a nonlinear constraint to an ML model.  The constraint
##' keeps the thresholds in order.  Constraints often slow model estimation,
##' however, keeping the thresholds in increasing order helps ensure the likelihood
##' function is well-defined.  If you're having problems with ordinal data, this is
##' one of the things to try.
##' @return
##' a new MxModel object with the constraints added
##' @seealso
##' \code{demo("omxConstrainMLThresholds")}
omxConstrainMLThresholds <- function(model, dist=.1) {
    expect <- model$expectation
    
	thresholdName <- expect$thresholds
    if(!is.na(thresholdName)) {
        threshData <- model$data
        if(suppressWarnings(is.na(model$data)) || 
            threshData$type != "raw") {
            stop(paste("Raw data must be in the model for",
                "omxConstrainThresholds to work."))
        }
        threshMat <- model[[thresholdName]]
        varNames <- colnames(threshMat)
        if(is.null(varNames)) {
            warning("No threshold names: defaulting to all data names")
            varNames <- names(threshData$observed)
        }
        for(i in seq_along(varNames)) {
            thisColumn <- threshData$observed[,varNames[i]]
            if(!is.ordered(thisColumn)) {
                stop(paste("Non-factor data found; could not constrain",
                        "column", i, "(", varNames[i], ")",
                        "of model", model$name))
            }
            nThresh <- nlevels(thisColumn)-1
            if(nThresh <= 1) {
                next;
            }
            conString <- paste0(thresholdName, "[1:",(nThresh-1), ",", i, "]", 
                                "<", thresholdName, "[2:", nThresh, ",", i, "]")
            conName <- paste0("ThresholdConstraint",i) # FIXME: May conflict!
            tCon <- eval(substitute(mxConstraint(theExpression,name = conName),
                list(theExpression = parse(text=conString)[[1]])))
            model <- mxModel(model, tCon)
        }
    }
    model
}
