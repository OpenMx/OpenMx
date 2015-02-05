omxConstrainMLThresholds <- function(model, dist=.1) {
    expect <- model$expectation
    
	thresholdName <- expect$thresholds
    if(!is.na(thresholdName)) {
        threshData <- model$data
        if(suppressWarnings(is.na(oneFactorJointModel$data)) || 
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
            nThresh <- nlevels(tDats)-1
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
