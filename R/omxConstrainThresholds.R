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
##' @examples
##' # a simple one factor ordinal model
##' require(OpenMx)
##' 
##' data(myFADataRaw)
##' 
##' oneFactorOrd <- myFADataRaw[,c("z1", "z2", "z3")]
##'
##' oneFactorOrd$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
##' oneFactorOrd$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
##' oneFactorOrd$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))
##' 
##'	oneFactorModel <- mxModel("Common Factor Model Path Specification", 
##'	type="RAM",
##'	mxData(
##'		observed=oneFactorOrd,
##'		type="raw"
##'	),
##'	manifestVars=c("z1","z2","z3"),
##'	latentVars="F1",
##'	# residual variances
##'	mxPath(
##'		from=c("z1","z2","z3"),
##'		arrows=2,
##'		free=FALSE,
##'		values=c(1,1,1),
##'		labels=c("e1","e2","e3")
##'	),
##'	# latent variance
##'	mxPath(
##'		from="F1",
##'		arrows=2,
##'		free=TRUE,
##'		values=1,
##'		labels ="varF1"
##'	),
##'	# factor loadings
##'	mxPath(
##'		from="F1",
##'		to=c("z1","z2","z3"),
##'		arrows=1,
##'		free=c(FALSE,TRUE,TRUE),
##'		values=c(1,1,1),
##'		labels=c("l1","l2","l3")
##'	),
##'	# means
##'	mxPath(
##'		from="one",
##'		to=c("z1","z2","z3","F1"),
##'		arrows=1,
##'		free=FALSE,
##'		values=0,
##'		labels=c("meanz1","meanz2","meanz3","meanF")
##'	),
##'	# thresholds
##'	mxThreshold(vars=c("z1", "z2", "z3"),
##'		nThresh=c(1,1,2),
##'		free=TRUE,
##'		values=c(-1, 0, -.5, 1.2)
##'		)
##') # close model
##' 
##' oneFactorCon <- omxConstrainMLThresholds(oneFactorModel)
##' #oneFactorResults <- mxRun(oneFactorCon)
##' #N.B. FAILS!
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
