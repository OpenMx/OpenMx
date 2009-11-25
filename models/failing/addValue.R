# SCRIPT: addValue.R
# Author: Timothy  Bates tim.bates@ed.ac.uk
# History:  Wed Nov 25 12:35:51 GMT 2009
# / modified omx functions to get get values onto paths
#  would be rolled into the existing functions when debugged
# Rhelp: http://www.statmethods.net
# OpenMx: http://www.openmx.virginia.com
##########################################
require(OpenMx); require(foreign); require(MASS); require(ggplot2); library(reshape);


# BEGIN Test model
  # Fig2.12 from John Loehlin edition 4
  # example pdf
  # http://www.subjectpool.com/data/2009/andrea/loehlin2.12.pdf
  require(OpenMx)
  names     = c("A","C","D")
  dimnames  = list(r=names, c=names)
  data      = matrix(c(1, .2,.24, .2, 1, .3, .24, .3, 1), nrow=3, byrow=TRUE, dimnames=list(r=names, c=names))
  manifests = row.names(data)
  latents   = c("B")

  model = mxModel("prob2-12", type="RAM", 
    manifestVars = manifests,
    latentVars   = latents,
    mxPath(from = "A", to="B", arrows=1), 
    mxPath(from = "B", to= c("C", "D"), arrows=1 ), 
    mxPath(from = c("A", "B", "C", "D"), arrows=2),   # residual variances
    mxData(data, type="cov",numObs=100)
  )

  fit = mxRun(model)
  myOmxGraphviz(fit , "test.dot"); system("open test.dot")
# END of test system

# Modified omx functions (prefixed with "my")

myOmxGraphviz = function (model, dotFilename = "") {
    if (missing(model) || !is(model, "MxModel")) {
        stop("The first argument is not an MxModel object")
    }
    if (!is(model, "MxRAMModel")) {
        stop(paste("The model", omxQuotes(model@name), "is not a 'RAM' type model"))
    }
    graph <- new("MxRAMGraph", model@manifestVars, model@latentVars)
    uniPaths    <- myMatrixToPaths(model[["A"]], 1)
    biPaths     <- myMatrixToPaths(model[["S"]], 2)
    meanPaths   <- meansToPaths(model[["M"]]) # todo: needs changin too??
    graph@paths <- c(graph@paths, uniPaths, biPaths, meanPaths)
    myWriteDotFile(model, graph, dotFilename)
}


myMatrixToPaths = function (mxMatrix, arrows = c(1, 2)) {
    values <- mxMatrix@values
    free   <- mxMatrix@free
    labels <- mxMatrix@labels
    select <- (values != 0) | (free) | (!is.na(labels))
    if (length(select) > 0) {
        rowFactors <- row(values, as.factor = TRUE)
        colFactors <- col(values, as.factor = TRUE)
        fromNames  <- as.character(colFactors[select])
        toNames    <- as.character(rowFactors[select])
        # Tue Nov 24 18:44:08 GMT 2009: (tb) modified to include theValues
        theValues  <- as.character(values[select])
        # return(theValues)
        if (length(fromNames) > 0 && length(toNames) > 0) {
          if (length(theValues) > 0){
            # Tue Nov 24 18:44:08 GMT 2009: (tb) modified to add theValues 
            return(mxPath(from = fromNames, to = toNames, arrows = arrows, values=theValues))
          } else {
            return(mxPath(from = fromNames, to = toNames, arrows = arrows, values=theValues))
          }
        }
    }
    return(list())
}


myWriteDotFile = function (model, graph, dotFilename) {
    outputLines <- list()
    graphName   <- paste("\"", model@name, "\"", sep = "")
    outputLines <- c(outputLines, paste("digraph", graphName, "{"))
    outputLines <- c(outputLines, paste("\t", "node [style=filled, fontname=\"Arial\", fontsize=16];"))
    if (length(graph@manifestVars) > 0) {
        for (i in 1:length(graph@manifestVars)) {
            outputLines <- c(outputLines, paste("\t", graph@manifestVars[[i]], 
                 "[shape=box, fillcolor=\"#a9fab1\", height=0.5, width=0.5];"))
        }
    }
    if (length(graph@latentVars) > 0) {
        for (i in 1:length(graph@latentVars)) {
            outputLines <- c(outputLines, paste("\t", graph@latentVars[[i]], "[shape=circle, fillcolor=\"#f4fd78\"];"))
        }
    }
    if (!is.null(model[["M"]])) {
        outputLines <- c(outputLines, paste("\t", "one", "[shape=triangle];"))
    }
    if (length(graph@paths) > 0) {
        for (i in 1:length(graph@paths)) {
            outputArrow <- list()
            path <- graph@paths[[i]]
            theLabel =  sub("^0.", ".", round(path$value,2)) # Tue Nov 24 18:22:00 GMT 2009: (tb) get values to use for path labels
            # Supressing the leading zero
            # TODO need to add parameter for rounding, and offset from line (just using a space below to force value off the line) ?
            outputArrow <- c(outputArrow, paste("\t", path$from, "->", path$to))
            if (path$arrows == 1) {
              # Tue Nov 24 18:22:00 GMT 2009: (tb) Added values to paths
              # TODO rewrite .dot to avoid cludge of adding a space to offset the label from the path
                outputArrow <- c(outputArrow, paste('[dir=forward, label =" ', theLabel, '"]', sep=""))
            } else if (path$arrows == 2) {
                if (path$from == path$to && path$from %in% graph@latentVars) {
                  # TODO: Would be nice to detect when a manifest var has a "to" connection with a 
                  #         latent var, and put the variance on top.
                  # TODO: Might lead to nicer layouting to leave the tailport info off... maybe a stick arrow
                  # TODO: Can graphviz let us say tailport=headport, and then the 
                  #         layouter figures the best place for the circle?
                  # TODO would be good for instruct graphviz to keep one-arrow lines straight
                  # Tue Nov 24 18:22:00 GMT 2009: (tb) Added values to paths
                  outputArrow <- c(outputArrow, paste('[dir=both, headport=n, tailport=n, label ="', theLabel, '"]', sep=""))
                }else if (path$from == path$to) {
                  # Tue Nov 24 18:22:00 GMT 2009: (tb) Added values to paths
                  outputArrow <- c(outputArrow, paste('[dir=both, headport=s, tailport=s, label ="', theLabel, '"]', sep=""))
                }else {
                  # Tue Nov 24 18:22:00 GMT 2009: (tb) Added values to paths
                  outputArrow <- c(outputArrow, paste("[dir=both, label =", theLabel, "]", sep=""))
                }
            }
            outputArrow <- c(outputArrow, ";")
            outputArrow <- paste(outputArrow, collapse = "")
            outputLines <- c(outputLines, outputArrow)
        }
    }
    outputLines <- c(outputLines, "}\n")
    outputLines <- paste(outputLines, collapse = "\n")
    cat(outputLines, file = dotFilename)
    return(invisible(outputLines))
}