#
#   Copyright 2007-2017 The OpenMx Project
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

mxCompare <- function(base, comparison, ..., all = FALSE,
		      boot=FALSE, replications=400, previousRun=NULL, checkHess=FALSE) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxCompare does not accept values for the '...' argument")
	}
	if (missing(base)) {
		stop("'base' argument be a MxModel object or list of MxModel objects")	
	}
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxCompare does not accept values for the '...' argument")
	}
	if (is.list(base)) {
		base <- unlist(base)
	} else {
		base <- list(base)
	}
	if(!all(sapply(base, is, "MxModel"))) {
		stop("The 'base' argument must consist of MxModel objects")
	}
	
	if (missing(comparison)) {
		# no comparison models, just make a dummy list to feed to showFitStatistics
		comparison <- list()
	} else {
		if (is.list(comparison)) {
			comparison <- unlist(comparison)
		} else {
			comparison <- list(comparison)
		}
		if(!all(sapply(comparison, is, "MxModel"))) {
			stop("The 'comparison' argument must consist of MxModel objects")
		}
	}
	if (missing(checkHess)) checkHess <- as.logical(NA)
	if (missing(boot) && (!missing(replications) || !missing(previousRun))) boot <- TRUE
	resultsTable <- showFitStatistics(base, comparison, all, boot, replications, previousRun, checkHess)
	return(resultsTable)
}

mxCompareMatrix <- function(models, statistic, ...,
			    boot=FALSE, replications=400, previousRun=NULL, checkHess=FALSE) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxCompareMatrix does not accept values for the '...' argument")
	}
	if (is.list(base)) {
		base <- unlist(base)
	} else {
		base <- list(base)
	}
	if(!all(sapply(base, is, "MxModel"))) {
		stop("The 'base' argument must consist of MxModel objects")
	}
	
	if (missing(checkHess)) checkHess <- as.logical(NA)
	if (missing(boot) && (!missing(replications) || !missing(previousRun))) boot <- TRUE

	resultsTable <- iterateNestedModels(models, boot, replications, previousRun, checkHess)

	if (missing(statistic)) stop(paste("Available statistics are:", omxQuotes(colnames(resultsTable))))

	unrecog <- !(statistic %in% c('raw', colnames(resultsTable)))
	if (any(unrecog)) {
		stop(paste("Statistic", omxQuotes(statistic[unrecog]), "is not available"))
	}
	if (statistic == 'raw') return(resultsTable)

	diagstat <- c('ep','minus2LL','df','AIC')
	offdiag <- c('diffLL','diffdf','p')

	modelNames <- sapply(models, function(m) m$name)
	mat <- matrix(NA, length(models), length(models),
		      dimnames=list(modelNames, modelNames))
	offdiagpick <- intersect(statistic, offdiag)
	if (length(offdiagpick)) {
		mat[,] <- resultsTable[,offdiagpick[1]]
	}
	diagpick <- intersect(statistic, diagstat)
	if (length(diagpick)) {
		diag(mat) <- resultsTable[!is.na(resultsTable$ep), diagpick[1]]
	}
	if (!is.null(attr(resultsTable, "bootData"))) {
		attr(mat, "bootData") <- attr(resultsTable, "bootData")
	}
	class(mat) <- 'result.mxCompareMatrix'
	mat
}

print.result.mxCompareMatrix <- function(x,...) print(x[,])

anova.MxModel <- function(object, ...) {
	args <- list(...)
	boot <- FALSE

	replications <- 400L
	ap <- match('replications', names(args))
	if (!is.na(ap)) {
		boot <- TRUE
		replications <- args[[ap]]
		args <- args[-ap]
	}
	
	previousRun <- NULL
	ap <- match('previousRun', names(args))
	if (!is.na(ap)) {
		boot <- TRUE
		previousRun <- args[[ap]]
		args <- args[-ap]
	}

	checkHess <- as.logical(NA)
	ap <- match('checkHess', names(args))
	if (!is.na(ap)) {
		checkHess <- args[[ap]]
		args <- args[-ap]
	}

	ap <- match('boot', names(args))
	if (!is.na(ap)) {
		boot <- args[[ap]]
		args <- args[-ap]
	}

	comparison <- args
	if (length(comparison) == 0) {
		stop(paste("Compare model",omxQuotes(object$name),
			   "with which other models?"))
	}
	if(!all(sapply(comparison, is, "MxModel"))) {
		stop("The '...' argument must consist of MxModel objects")
	}
	models <- c(list(object), comparison)
	summaries <- lapply(models, summary)
	stats <- data.frame(m2ll=unlist(Map(function(x) x$Minus2LogLikelihood, summaries)),
			    df  =unlist(Map(function(x) x$degreesOfFreedom, summaries)))
	mo <- with(stats, order(df, m2ll))
	models <- models[mo]
	stats <- stats[mo,]

	bootData <- NULL
	if (boot) {
		bootTodo <- list()
		for (rx in 2:length(models)) {
			if (stats[rx-1,'m2ll'] >= stats[rx,'m2ll'] ||
			    stats[rx-1,'df'] >= stats[rx,'df']) next
			nChar <- as.character(rx)
			bootTodo[[ nChar ]] <- rx-1
		}
		bootData <- setupBootData(models, models, bootTodo, replications, previousRun)
		bootData <- fillBootData(models, models, bootTodo, bootData, checkHess)
	}

	result <- list()

	result[[ length(result) + 1L ]] <- 
		collectBaseStatistics(newEmptyCompareRow(), models[[1]])

	for (i in 2:length(models)) {
		if (stats[i-1,'m2ll'] > stats[i,'m2ll']) {
			result[[ length(result) + 1L ]] <-
				collectBaseStatistics(newEmptyCompareRow(), models[[i]])
			next
		}
		boot1 <- extractLRTBootstrapPair(bootData, i, i-1)
		result[[ length(result) + 1L ]] <-
			collectStatistics(newEmptyCompareRow(), models[[i-1]], models[[i]], boot1)
	}

	ret <- do.call(rbind, result)
	if (boot) {
		attr(ret, "bootData") <- bootData
	}
	ret
}

newEmptyCompareRow <- function() {
	data.frame(stringsAsFactors = FALSE,
		   base=as.character(NA),
		   comparison=as.character(NA),
		   ep=as.numeric(NA),
		   minus2LL=as.numeric(NA),
		   df=as.numeric(NA),
		   AIC=as.numeric(NA),
		   diffLL=as.numeric(NA),
		   diffdf=as.numeric(NA), 
		   p=as.numeric(NA))
}

iterateNestedModels <- function(models, boot, replications, previousRun, checkHess) {
	result <- list()
	
	summaries <- lapply(models, summary)

	bd <- NULL
	if (boot) {
		bootTodo <- list()
		for (rx in seq_along(models)) {
			for (cx in seq_along(models)) {
				if (cx == rx ||
				    summaries[[rx]]$degreesOfFreedom > summaries[[cx]]$degreesOfFreedom) next
				cxChar <- as.character(cx)
				bootTodo[[cxChar]] <- c(bootTodo[[cxChar]], rx)
			}
		}

		bd <- setupBootData(models, models, bootTodo, replications, previousRun)
		bd <- fillBootData(models, models, bootTodo, bd, checkHess)
	}

	for (rx in seq_along(models)) {
		for (cx in seq_along(models)) {
			if (cx == rx) {
				r1 <- collectBaseStatistics(newEmptyCompareRow(), models[[cx]])
				result[[length(result) + 1]] <- r1
				next
			}
			rmod <- models[[rx]]
			cmod <- models[[cx]]
			if (summaries[[rx]]$degreesOfFreedom > summaries[[cx]]$degreesOfFreedom) {
				result[[length(result) + 1]] <- newEmptyCompareRow()
			} else {
				boot1 <- extractLRTBootstrapPair(bd, cx, rx)
				result[[length(result) + 1]] <-
					collectStatistics1(newEmptyCompareRow(), models[[rx]], models[[cx]], boot1)
			}
		}
	}

	ret <- do.call(rbind, result)
	if (boot) {
		attr(ret, "bootData") <- bd
	}
	ret
}

loadDataIntoModel <- function(model, dataList) {
  for (modelName in names(dataList)) {
    if (modelName == model$name) {
      model@data@observed <- dataList[[modelName]]
    } else {
      model[[modelName]]@data@observed <- dataList[[modelName]]
    }
  }
  model
}

setupBootData <- function(nullHyp, comparison, todo,
			  replications, previousRun) {
  
  # pre-check data compatibility of nullHyp and comparison? TODO

	if (!is.null(previousRun)) previousRun <- attr(previousRun, "bootData")

  seedVec <- as.integer(runif(replications, min = -2e9, max=2e9))
  
  bootData <- list()
  for (i in names(todo)) {
    nullHypData <- data.frame(seed=seedVec, fit=NA,
			      statusCode=as.statusCode(NA))
    for (par in names(coef( nullHyp[[ as.integer(i) ]] ))) nullHypData[[par]] <- NA
    prevData <- previousRun[[ i ]]
    if (!is.null(prevData)) {
	    if (ncol(prevData) != ncol(nullHypData) || any(colnames(prevData) != colnames(nullHypData))) {
		    warning("Data from previousRun does not match current models (ignored)")
	    } else {
		    rows <- min(nrow(prevData), replications)
		    nullHypData[1:rows,] <- prevData[1:rows,]
	    }
    }
    bootData[[ i ]] <- nullHypData

    for (j in todo[[i]]) {
      cmpData <- data.frame(seed=seedVec, fit=NA,
                            statusCode=as.statusCode(NA))
      for (par in names(coef( comparison[[j]] ))) cmpData[[par]] <- NA
      key <- paste(i,j,sep=":")
      prevData <- previousRun[[ key ]]
      if (!is.null(prevData)) {
	      if (ncol(prevData) != ncol(cmpData) || any(colnames(prevData) != colnames(cmpData))) {
		      warning("Data from previousRun does not match current models (ignored)")
	      } else {
		      rows <- min(nrow(prevData), replications)
		      cmpData[1:rows,] <- prevData[1:rows,]
	      }
      }
      bootData[[ key ]] <- cmpData
    }
  }
  bootData
}

fillBootData <- function(nullHyp, comparison, todo, bootData, checkHess) {
  replications <- nrow(bootData[[1]])
  numThreads <- imxGetNumThreads()
  if (numThreads < 2L) {
	  workPlan <- factor(rep(1,replications))
  } else {
	  workPlan <- cut(1:replications, numThreads)
  }
  slices <- omxLapply(levels(workPlan), function(lev) {
	  for (repl in which(workPlan == lev)) {
		  for (i in names(todo)) {
			  nullHypData <- bootData[[ i ]]
			  if (!is.na(nullHypData[repl, 'fit'])) next

			  set.seed(nullHypData[repl, 'seed'])
			  nullModel <- nullHyp[[ as.integer(i) ]]
			  simData <- try(mxGenerateData(nullModel, returnModel=FALSE), silent=TRUE)
			  if (is(simData, "try-error")) {
				  stop(paste("Cannot bootstrap null model", omxQuotes(nullModel$name),
					     "; Does this model contain raw data?"), call.=FALSE)
			  }
			  if (is(simData, "data.frame")) {
				  simData <- list(simData)
				  names(simData) <- nullModel$name
			  }
			  null1 <- loadDataIntoModel(nullModel, simData)
			  null1 <- ProcessCheckHess(null1, checkHess)
			  null1 <- mxRun(null1, silent=TRUE, suppressWarnings = TRUE)
			  nullHypData[repl, 'fit'] <- null1$output$fit
			  nullHypData[repl, names(coef(null1))] <- coef(null1)
			  nullHypData[repl, 'statusCode'] <- as.statusCode(null1$output$status$code)
			  
			  bootData[[ i ]] <- nullHypData

			  topDataIndex <- match(null1$name, names(simData))
			  
			  for (j in todo[[i]]) {
				  cmpData <- bootData[[ paste(i,j,sep=":") ]]
				  if (!is.na(cmpData[repl, 'fit'])) next

				  cmpModel <- comparison[[j]]
				  names(simData)[topDataIndex] <- cmpModel$name
				  cmp1 <- loadDataIntoModel(cmpModel, simData)
				  cmp1 <- ProcessCheckHess(cmp1, checkHess)
				  cmp1 <- mxRun(cmp1, silent=TRUE, suppressWarnings = TRUE)

				  cmpData[repl, 'fit'] <- cmp1$output$fit
				  cmpData[repl, names(coef(cmp1))] <- coef(cmp1)
				  cmpData[repl, 'statusCode'] <- as.statusCode(cmp1$output$status$code)
				  
				  bootData[[ paste(i,j,sep=":") ]] <- cmpData
			  }
		  }
	  }
	  bootData
  })
  names(slices) <- levels(workPlan)
  for (lev in levels(workPlan)) {
	  slice <- slices[[lev]]
	  for (boot1 in names(slice)) {
		  bootData[[boot1]][workPlan == lev,] <- slice[[boot1]][workPlan == lev,]
	  }
  }
  bootData
}

extractLRTBootstrapPair <- function(bootData, i, j) {
	if (is.null(bootData)) return(NULL)
	boot1 <- bootData[[ as.character(i) ]]
	boot2 <- bootData[[ paste(i,j,sep=":") ]]
	if (!is.null(boot1) && !is.null(boot2)) return(list(boot1,boot2))
	NULL
}

showFitStatistics <- function(base, compare, all, boot, replications, previousRun, checkHess)  {
	bootData <- NULL
	if (boot) {
		bootTodo <- list()
		if (all) {
			for (rx in seq_along(base)) {
				for (cx in seq_along(compare)) {
					cxChar <- as.character(cx)
					bootTodo[[cxChar]] <- c(bootTodo[[cxChar]], rx)
				}
			}
		} else {
			maxLength <- max(length(base), length(compare))
			for (i in 1:maxLength) {
				baseIndex <- (i - 1L)%%length(base) + 1L
				compareIndex <- (i - 1L)%%length(compare) + 1L
				ciChar <- as.character(compareIndex)
				bootTodo[[ciChar]] <- c(bootTodo[[ciChar]], baseIndex)
			}
		}

		bootData <- setupBootData(compare, base, bootTodo, replications, previousRun)
		bootData <- fillBootData(compare, base, bootTodo, bootData, checkHess)
	}

    statistics <- list()
    if (all) {
        for (i in seq_along(base)) {
            nextBase <- base[[i]]
            statistics[[length(statistics) + 1]] <-
		    collectBaseStatistics(newEmptyCompareRow(), nextBase)
            for (j in seq_along(compare)) {
                nextCompare <- compare[[j]]
		boot1 <- extractLRTBootstrapPair(bootData, i, j)
                statistics[[length(statistics) + 1]] <-
			collectStatistics(newEmptyCompareRow(), nextBase, nextCompare, boot1)
            }
        }
    }
    else {
		if(length(compare)==0){
	        for (i in seq_along(base)) {
			statistics[[length(statistics) + 1]] <-
				collectBaseStatistics(newEmptyCompareRow(), base[[i]])
	        }
		} else {
	        maxLength <- max(length(base), length(compare))
	        previousBaseIndex <- -1
	        for (i in 1:maxLength) {
	            nextBaseIndex <- (i - 1)%%length(base) + 1
	            nextCompareIndex <- (i - 1)%%length(compare) + 1
	            nextBase <- base[[nextBaseIndex]]
	            nextCompare <- compare[[nextCompareIndex]]
	            if (previousBaseIndex != nextBaseIndex) {
	                statistics[[length(statistics) + 1]] <-
				collectBaseStatistics(newEmptyCompareRow(), nextBase)
	            }
		    boot1 <- extractLRTBootstrapPair(bootData, nextBaseIndex, nextCompareIndex)
	            statistics[[length(statistics) + 1]] <-
			    collectStatistics(newEmptyCompareRow(), nextBase, nextCompare, boot1)
	            previousBaseIndex <- nextBaseIndex
	        }
		}		
    }
	statistics <- do.call(rbind, statistics)
	if (boot) {
		attr(statistics, "bootData") <- bootData
	}
    return(statistics)
}

collectBaseStatistics <- function(row, ref) {
	refSummary <- summary(ref)
	row[,'base'] <- refSummary$modelName
	row[,c('ep','minus2LL','df','AIC')] <-
		c(refSummary$estimatedParameters,
		  refSummary$Minus2LogLikelihood,
		  refSummary$degreesOfFreedom,
		  refSummary$AIC.Mx)
	row
}

collectStatistics <- function(otherStats, ref, other, bootPair) {
	otherStats <- collectBaseStatistics(otherStats, other)
	collectStatistics1(otherStats, ref, other, bootPair)
}

collectStatistics1 <- function(otherStats, ref, other, bootPair) {
	refSummary <- summary(ref)
	otherSummary <- summary(other)
	otherStats[,c('base','comparison')] <-
		c(refSummary$modelName,
		  otherSummary$modelName)
	otherStats[,c('diffLL','diffdf')] <-
		c(otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood,
		  otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom)

	if (any(otherStats[['diffdf']] < 0)) {
		msg <- paste("Model", omxQuotes(refSummary$modelName), "has more degrees of freedom than",
			     otherSummary$modelName, "which means that the models need to be",
			     "compared in the opposite order")
		warning(msg)
	} else {
		if (is.null(bootPair)) {
			otherStats[['p']] <- pchisq(otherStats[['diffLL']], otherStats[['diffdf']], lower.tail=FALSE)
		} else {
			baseData <- bootPair[[1]] # null hypothesis
			cmpData <- bootPair[[2]]
			
			if (any(baseData[,'seed'] != cmpData[,'seed'])) {
				stop("Some seeds do not match")
			}
			
			# If optimization goes wrong, the null model could fit better than
			# the alternative model. We exclude these replications.
			mask <- (baseData[,'statusCode'] %in% mxOption(other, "Status OK") &
				 cmpData[,'statusCode'] %in% mxOption(ref, "Status OK") &
				 baseData[,'fit'] - cmpData[,'fit'] > 0)
			if (sum(mask) < 3) {
				stop(paste("Fewer than 3 replications are available."))
			}
			if (sum(mask) < .95*length(mask)) {
				pct <- round(100*sum(mask) / length(mask))
				warning(paste0("Only ",pct,"% of the bootstrap replications ",
					       "converged acceptably. Accuracy is much less than the ", nrow(raw),
					       " replications requested"), call.=FALSE)
			}

			diffLL <- baseData[mask,'fit'] - cmpData[mask,'fit']
			otherStats[['p']] <- sum(diffLL > otherStats[['diffLL']]) / sum(mask)
		}
		otherStats[['p']][ otherStats[['diffdf']] == 0 ] <- NA
	}
	return(otherStats)
}

mxParametricBootstrap <- function(nullModel, labels,
                                  alternative=c("two.sided", "greater", "less"),
				  ...,
                                  alpha=0.05,
                                  correction=p.adjust.methods,
                                  previousRun=NULL, replications=400,
                                  checkHess=FALSE,
				  signif.stars = getOption("show.signif.stars")) {
	garbageArguments <- list(...)
	if (length(garbageArguments) > 0) {
		stop("mxParametricBootstrap does not accept values for the '...' argument")
	}
	if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
	if (missing(checkHess)) checkHess <- as.logical(NA)
	nullModel <- ProcessCheckHess(nullModel, checkHess)
  alternative <- match.arg(alternative)
  correction <- match.arg(correction)
  seedVec <- as.integer(runif(replications, min = -2e9, max=2e9))
  bootData <- data.frame(seed=seedVec, fit=NA,
                         statusCode=as.statusCode(NA))
  for (par in c(names(coef(nullModel)), labels)) bootData[[par]] <- NA
  
	if (!is.null(previousRun)) {
		previousRun <- attr(previousRun, "bootData")
		if (ncol(previousRun) != ncol(bootData) || any(colnames(previousRun) != colnames(bootData))) {
			warning("Data from previousRun does not match current model (ignored)")
		} else {
			rows <- min(nrow(previousRun), replications)
			bootData[1:rows,] <- previousRun[1:rows,]
		}
	}

  for (repl in 1:replications) {
    set.seed(bootData[repl, 'seed'])
    if (!is.na(bootData[repl, 'fit'])) next
    base1 <- mxGenerateData(nullModel, returnModel=TRUE)
    base1 <- omxSetParameters(base1, labels=labels, free=TRUE)
    base1 <- mxRun(base1, silent=TRUE, suppressWarnings = TRUE)
    bootData[repl, 'fit'] <- base1$output$fit
    bootData[repl, names(coef(base1))] <- coef(base1)
    bootData[repl, 'statusCode'] <- as.statusCode(base1$output$status$code)
  }
  
  nullParam <- omxGetParameters(nullModel, free=NA)[labels]
  model <- omxSetParameters(nullModel, labels=labels, free=TRUE)
  model <- mxRun(model, silent = TRUE, suppressWarnings = TRUE)
  
  mask <- bootData[,'statusCode'] %in% mxOption(nullModel, "Status OK")
  if (sum(mask) < 3) {
    stop(paste("Fewer than 3 replications are available."))
  }
  if (sum(mask) < .95*length(mask)) {
    pct <- round(100*sum(mask) / length(mask))
    warning(paste0("Only ",pct,"% of the bootstrap replications ",
                   "converged acceptably. Accuracy is much less than the ",
                   replications,
                   " replications requested"), call.=FALSE)
  }
  
  est <- coef(model)[labels]
  pval <- rep(NA, length(labels))
  for (lx in seq_along(labels)) {
    if (alternative == "two.sided") {
      pval[lx] <- sum(abs(bootData[mask,labels[lx]]) > abs(est[labels[lx]]))
    } else if (alternative == "greater") {
      pval[lx] <- sum(bootData[mask,labels[lx]] > est[labels[lx]])
    } else {  # less
      pval[lx] <- sum(bootData[mask,labels[lx]] < est[labels[lx]])
    }
  }
  pval <- pval / sum(mask)
  
  pzero <- pval == 0
  pval[pzero] <- 1/sum(mask)
  pval <- p.adjust(pval, method=correction)

  ret <- data.frame(est=est,
             null=nullParam,
             p=pval)
  if (signif.stars) {
    sig <- rep('', length(labels))
    sig[pval < alpha] <- '*'
    sig[pval < alpha/5] <- '**'
    sig[pval < alpha/50] <- '***'
    ret$sig <- sig
  }
  ret$note <- ''
  ret[pzero,'note'] <- paste0('< 1/',sum(mask))
  attr(ret,'bootData') <- bootData
  ret
}
