#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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
	prohibitDotdotdot(list(...))
	if (missing(base)) {
		stop("'base' argument be a MxModel object or list of MxModel objects")	
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

mxCompareMatrix <- function(models,
			    diag=c('minus2LL','ep','df','AIC'),
			    stat=c('p', 'diffLL','diffdf'), ...,
			    boot=FALSE, replications=400, previousRun=NULL, checkHess=FALSE,
			    wholeTable=FALSE) {
	prohibitDotdotdot(list(...))
	
	if (missing(checkHess)) checkHess <- as.logical(NA)
	if (missing(boot) && (!missing(replications) || !missing(previousRun))) boot <- TRUE

	diagpick <- match.arg(diag)
	offdiagpick <- match.arg(stat)

	resultsTable <- iterateNestedModels(models, boot, replications, previousRun, checkHess)

	if (wholeTable) return(resultsTable)

	modelNames <- sapply(models, function(m) m$name)
	mat <- matrix(NA, length(models), length(models),
		      dimnames=list(modelNames, modelNames))
	mat[,] <- resultsTable[,offdiagpick[1]]
	diag(mat) <- resultsTable[!is.na(resultsTable$ep), diagpick[1]]
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

assertIsRawData <- function(model) {
  type <- model$data$type
  if (type == 'raw') return()
  stop(paste("Model", omxQuotes(model$name), "contains", omxQuotes(type),
	     "data. Only type='raw' data is supported by",
	     "mxPowerSearch(..., method='empirical')"))
}

loadDataIntoModel <- function(model, dataList, assertRaw=FALSE) {
  for (modelName in names(dataList)) {
    dataobj <- mxData(dataList[[modelName]], type='raw')
    if (modelName == model$name) {
      if (assertRaw) assertIsRawData(model)
      model <- mxModel(model, dataobj)
    } else {
      if (assertRaw) assertIsRawData(model[[modelName]])
      model <- mxModel(model, mxModel(model[[modelName]], dataobj))
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
			  simData <- try(mxGenerateData(nullModel, returnModel=FALSE))
			  if (is(simData, "try-error")) {
				  stop(paste("Cannot bootstrap null model", omxQuotes(nullModel$name)), call.=FALSE)
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
	assertModelRunAndFresh(ref)
	assertModelRunAndFresh(other)
	refSummary <- summary(ref)
	otherSummary <- summary(other)
	
	#Check for validity of the comparison ###
	rfu <- ref$output$fitUnits #<--NULL if output slot is empty.
	if(!length(rfu)){
		warning(paste("MxModel '",ref$name,"' has no 'fitUnits' element in its output slot; has it been run?",sep=""))
	}
	ofu <- other$output$fitUnits #<--NULL if output slot is empty.
	if(!length(ofu)){
		warning(paste("MxModel '",other$name,"' has no 'fitUnits' element in its output slot; has it been run?",sep=""))
	}
	#Only stop if there's a definite mismatch in fit units:
	if(length(rfu) && length(ofu) && rfu!=ofu){
		stop(paste("MxModel '",ref$name,"' has '",rfu,"' fit units, but MxModel '",other$name,"' has '",ofu,"' fit units",sep=""))
	}
	#Even though the fit units match, the restricted ML and ordinary ML fit values can't be validly compared:
	if( is(ref$fitfunction,"MxFitFunctionGREML")!=is(other$fitfunction,"MxFitFunctionGREML") ){
		stop(paste("MxModel '",ref$name,"' has a fitfunction of class '",class(ref$fitfunction),"', but MxModel '",other$name,"' has a fitfunction of class '",class(other$fitfunction),"'",sep=""))
	}
	rgfe <- refSummary$GREMLfixeff #<--NULL unless model uses a GREML expectation and has been run
	if(length(rgfe)){rgfe <- paste(rgfe$name,collapse=",")}
	ogfe <- otherSummary$GREMLfixeff #<--NULL unless model uses a GREML expectation and has been run
	if(length(ogfe)){ogfe <- paste(ogfe$name,collapse=",")}
	if( length(rgfe)!=length(ogfe) || (length(rgfe) && length(ogfe) && rgfe!=ogfe) ){
		#This is a warning, not an error, because it's possible that the user is indeed using the same covariates in both models, but with
		#different column names.  (If one of the models hasn't been run yet, GREMLFixEffList() will return NULL, but the fit value will 
		#be NA, so the output for the comparison won't even look valid):
		warning(paste("the names of the fixed effects in MxModels '",ref$name,"' and '",other$name,"' do not match; comparison of REML fit values is only valid for models that use the same covariates",sep=""))
	}
	#End validity checks
	
	otherStats[,c('base','comparison')] <-
		c(refSummary$modelName,
		  otherSummary$modelName)
	otherStats[,c('diffLL','diffdf')] <-
		c(otherSummary$Minus2LogLikelihood - refSummary$Minus2LogLikelihood,
		  otherSummary$degreesOfFreedom - refSummary$degreesOfFreedom)

  diffdf <- otherStats[['diffdf']]
  diffdf <- diffdf[!is.na(diffdf)]
	if (length(diffdf) && any(diffdf < 0)) {
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
	prohibitDotdotdot(list(...))
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

totalSampleSize <- function(model, default=100L) {
  sizes <- sapply(extractData(model), function(mxd) {
    if (mxd@type == 'raw') {
      nrow(mxd@observed)
    } else {
      mxd@numObs
    }
  })
  if (length(sizes) == 0) sizes <- default
  sum(sizes)
}

fitPowerModel <- function(rx, result) {
  # rx is which(is.na(result$reject))[1] - 1L
  result <- result[!is.na(result$x),]
  algRle <- rle(result$alg)
  if (algRle$values[1] == 'init') {
    # skip outliers from early probe sequence
    newFirst <- max(1, algRle$lengths[1] - 3L)
    result <- result[newFirst:nrow(result),]
  }
  m1 <- suppressWarnings(glm(reject ~ x, data = result, family = binomial))
  if (all(!is.na(coef(m1))) && summary(m1)$coefficients['x','Pr(>|z|)'] < .25) {
    curX <- as.numeric(((rx %% 3) - coef(m1)[1]) / coef(m1)[2])
    alg <- '2p'
  } else {
    m2 <- glm(reject ~ 1, data = result, family = binomial)
    from <- max(nrow(result)-9, 1)
    curX <- mean(result$x[from:nrow(result)]) * ifelse(coef(m2)[1] < 0, 1.1, 0.9)
    alg <- '1p'
  }
    # Can only consider one-sided hypotheses
  if (sign(curX) != sign(result[1,'x'])) curX <- 0
  list(curX=curX, m1=m1, alg=alg)
}

validateSigLevel <- function(sl) {
  if (length(sl) > 1) {
    stop(paste("Can only evaluate one sig.level at a time.",
	       "To evaluate power across a range of sig.levels",
	       "you need to call mxPower in a loop and accumulate",
	       "estimates that way"))
  } else if (length(sl) == 0) {
    stop("At what sig.level?")
  } else if (sl <= 0 || sl >= 1) {
    stop("sig.level must be between 0 and 1")
  }
}

mxPowerSearch <- function(trueModel, falseModel, n=NULL, sig.level=0.05, ...,
	probes = 300L, previousRun = NULL,
	gdFun = mxGenerateData,
	method = c('empirical', 'ncp'),
	grid = NULL,
	statistic = c('LRT','AIC','BIC'),
	OK = mxOption(trueModel, "Status OK"), checkHess=FALSE,
	silent = !interactive())
{
	# TODO: add plot=TRUE? or return S3 object that responds to plot(obj)?
    garbageArguments <- list(...)
    if (length(garbageArguments) > 0) {
		 message("Invalid inputs to mxPowerSearch:")
		 print(garbageArguments)
       stop("mxPowerSearch does not accept values for the '...' argument\n")
    }
  origSampleSize <- totalSampleSize(trueModel)
  validateSigLevel(sig.level)
    method <- match.arg(method)
    statistic <- match.arg(statistic)
    if (method == 'ncp') {
    if (!is.null(n)) stop(paste("method='ncp' does not work for fixed n =", n))
    if (statistic != 'LRT') stop(paste("method='ncp' does not work for statistic =", statistic))
    warnModelCreatedByOldVersion(trueModel)
    warnModelCreatedByOldVersion(falseModel)
    if (!trueModel@.wasRun || trueModel@.modifiedSinceRun) {
		warning("Polite but strong warning: You haven't re-run trueModel since modifying it.\n",
		"You most likely want to mxRun trueModel, or your apparent power may be illusory.\n",
		"(But perhaps you're an expert wanting to do this - if so we're glad you chose OpenMx for your advanced studies :-) ).")
	    trueModel <- mxRun(mxModel(trueModel, mxComputeOnce('fitfunction','fit')))
    }
    if (!falseModel@.wasRun || falseModel@.modifiedSinceRun) {
		warning("Polite but strong warning: You haven't re-run falseModel since modifying it.\n",
		"You most likely want to mxRun falseModel, or your apparent power may be illusory.\n",
		"(But perhaps you're an expert wanting to do this - if so we're glad you chose OpenMx for your advanced studies :-) ).")
	    falseModel <- mxRun(mxModel(falseModel, mxComputeOnce('fitfunction','fit')))
    }
    if (trueModel$output[['fitUnits']] != '-2lnL') {
	    stop(paste(trueModel$name, "measured in terms of", trueModel$output[['fitUnits']],
		       "instead of -2lnL"))
    }
    if (falseModel$output[['fitUnits']] != '-2lnL') {
	    stop(paste(falseModel$name, "measured in terms of", falseModel$output[['fitUnits']],
		       "instead of -2lnL"))
    }
    avgNcp <- falseModel$output$Minus2LogLikelihood - trueModel$output$Minus2LogLikelihood
    if (avgNcp < 0) stop("falseModel fit better than trueModel?")
    if (avgNcp == 0.0) stop("falseModel and trueModel are identical?")
    # is diffdf>1 ever a good approx? TODO
    diffdf <- summary(falseModel)[['degreesOfFreedom']] - summary(trueModel)[['degreesOfFreedom']]
    
    if (is.null(grid)) {
      width <- 2.75/avgNcp
      center <- 1.15*qchisq(1 - sig.level, diffdf)/avgNcp
      grid <- round(seq(center-1*width,
			center+4*width, length.out = 20) * origSampleSize)
    }
    out <- data.frame(x=grid)
    out$power <- 1 - suppressWarnings(pchisq(qchisq(1 - sig.level, diffdf), diffdf,
					     avgNcp * out$x / origSampleSize))
    if (any(is.na(out$power))) {
      stop(paste("Sorry, unable to estimate power for sig.level=",sig.level,
		 "with method='ncp'. Try method='empirical'"))
    }
    out$lower <- NA
    out$upper <- NA
    colnames(out)[1] <- 'N'
    return(out)
  }
  
  interest <- setdiff(names(coef(trueModel)), names(coef(falseModel)))
  if (!is.null(n) && length(interest) != 1) {
    stop(paste("Specify only 1 parameter (not", omxQuotes(interest),
               ") to search a parameter:power relationship"))
  }
  xLabel <- ''
  if (is.null(n)) {
    xLabel <- 'N'
    if (!silent) message(paste("Search n:power relationship for", omxQuotes(interest)))
  } else {
    xLabel <- interest
    if (!silent) message(paste0("Search ",interest,":power relationship for n=", n))
  }

  result <- data.frame(seed=as.integer(runif(probes, min = -2e9, max=2e9)),
		       reject=NA, x=NA, alg=NA, mseTrue=NA,
		       statusTrue=as.statusCode(NA), statusFalse=as.statusCode(NA))
  if (is.null(n)) {
    nullInterestValue <- 0
    curX <- 1.0
  } else {
    par <- omxGetParameters(falseModel, free=FALSE, labels=interest)
    if (!(interest %in% names(par))) {
      stop(paste("Cannot find", omxQuotes(interest),
                 "in falseModel. Please label it in both models and try again"))
    }
    nullInterestValue <- par[interest]
    curX <- (coef(trueModel)[interest] - nullInterestValue)/2
    if (curX == 0.0) curX <- .1
  }
  m1 <- NULL
  nextTrial <- 1L
  alg <- 'init'

  if (!is.null(previousRun)) {
      prevArgs <- attr(previousRun, "arguments")
      oldProbes <- attr(previousRun, 'probes')
      if (is.null(prevArgs$n) != is.null(n)) {
          warning("previousRun references a different kind of search (ignored)")
      } else if (!is.null(n) && prevArgs$n != n) {
          warning("previousRun searched a different sample size (ignored)")
      } else if (prevArgs$statistic != statistic) {
          warning("previousRun used a different statistic (ignored)")
      } else if (statistic == 'LRT' && prevArgs$sig.level != sig.level) {
          warning("previousRun used a different sig.level (ignored)")
      } else if (is.null(oldProbes)) {
        warning("previousRun did not contain old probes (ignored)")
      } else if (!all(colnames(oldProbes) == colnames(result))) {
        warning("previousRun old probes in wrong format (ignored)")
      } else {
        toCopy <- min(probes, nrow(oldProbes))
        result[1:toCopy,] <- oldProbes[1:toCopy,]
	okResult <- result[result$statusTrue %in% OK & result$statusFalse %in% OK,]
        nextTrial <- which(is.na(result$reject))[1]
	pm <- fitPowerModel(ifelse(!is.na(nextTrial), nextTrial-1L, nrow(result)), okResult)
	m1 <- pm$m1
	curX <- pm$curX
      }
  }

    trueModel <- ProcessCheckHess(trueModel, checkHess)
    falseModel <- ProcessCheckHess(falseModel, checkHess)

  prevProgressLen <- 0L
  if (!is.na(nextTrial)) for (rx in nextTrial:probes) {
    set.seed(result[rx,'seed'])
    info <- paste0(xLabel,"[", rx,"] fitting model '", alg, "' value ",
		   ifelse(is.null(n), round(origSampleSize * (nullInterestValue + curX)),
			  nullInterestValue + curX))
    if (!silent) imxReportProgress(info, prevProgressLen)
    prevProgressLen <- nchar(info)
    if (!is.null(n)) {
      trueModel <- omxSetParameters(trueModel, labels=interest,
                                    values = nullInterestValue + curX)
    }
    simData <- try(gdFun(trueModel, returnModel=FALSE,
                         nrowsProportion=ifelse(is.null(n), curX, n/origSampleSize)))
    if (is(simData, "try-error")) {
      stop(paste("Cannot generate data with trueModel",
                 omxQuotes(trueModel$name)), call.=FALSE)
    }
    if (is(simData, "data.frame")) {
      simData <- list(simData)
      names(simData) <- trueModel$name
    }
    
    true1  <- loadDataIntoModel(trueModel,  simData, assertRaw=rx==1)
    true1  <- mxRun(true1,  silent=TRUE, suppressWarnings = TRUE)
    # complain about parameters at box constraints TODO
    
    topDataIndex <- match(trueModel$name, names(simData))
    names(simData)[topDataIndex] <- falseModel$name
    false1 <- loadDataIntoModel(falseModel, simData)
    false1 <- mxRun(false1, silent=TRUE, suppressWarnings = TRUE)
    
    if (statistic == 'LRT') {
	    cmp1 <- mxCompare(true1, false1)
	    pval <- cmp1[2,'p']
	    if (is.na(pval)) stop("falseModel fit better than trueModel?")
	    result[rx, 'reject'] <- pval < sig.level
    } else if (statistic == 'AIC') {
	    rej <- summary(true1)[['informationCriteria']]['AIC:','df'] <
		    summary(false1)[['informationCriteria']]['AIC:','df']
	    result[rx, 'reject'] <- rej
    } else if (statistic == 'BIC') {
	    rej <- summary(true1)[['informationCriteria']]['BIC:','df'] <
		    summary(false1)[['informationCriteria']]['BIC:','df']
	    result[rx, 'reject'] <- rej
    } else { stop(statistic) }
    result[rx, 'x'] <- curX
    result[rx, 'alg'] <- alg
    result[rx, 'mseTrue'] <- sum((coef(true1) - coef(trueModel))^2)
    result[rx, 'statusTrue'] <- as.statusCode(true1$output$status$code)
    result[rx, 'statusFalse'] <- as.statusCode(false1$output$status$code)

    okResult <- result[result$statusTrue %in% OK & result$statusFalse %in% OK,]
    rej <- table(okResult$reject)
    if (dim(rej) == 1 || any(rej < 2)) {
      if (names(sort(rej, decreasing=TRUE))[1] == "TRUE") {
        curX <- curX / 2
      } else {
        curX <- curX * 2
      }
      alg <- 'init'
    } else {
      pm <- fitPowerModel(rx, okResult)
      m1 <- pm$m1
      curX <- pm$curX
      alg <- pm$alg
    }
  }
  if (!silent) imxReportProgress('', prevProgressLen)
  if (is.null(m1)) stop("Logistic model failed to converge")
  if (is.null(grid)) {
    width <- 1/coef(m1)[2]
    center <- -(coef(m1)[1] / coef(m1)[2])
    grid <- seq(max(center-1*width,0), center+4*width, length.out = 20)
    if (is.null(n)) grid <- round(grid * origSampleSize)
  }
  out <- data.frame(x=grid)
#  out$p <- plogis(out$N, center, 1/coef(m1)[2])
  if (is.null(n)) {
    pr <- predict(m1, newdata=out / origSampleSize, type="link", se.fit=TRUE)
  } else {
    pr <- predict(m1, newdata=out, type="link", se.fit=TRUE)
  }
  out$power <- plogis(pr$fit)
  out$lower <- plogis(pr$fit - 2*pr$se.fit)
  out$upper <- plogis(pr$fit + 2*pr$se.fit)
  attr(out, "probes") <- result
  attr(out, "arguments") <- list(n=n, sig.level=sig.level, statistic=statistic)
  attr(out, "model") <- m1
  out$x <- out$x + nullInterestValue
  colnames(out)[1] <- xLabel
  out
}

mxPower <- function(trueModel, falseModel, n=NULL, sig.level=0.05, power=0.8, ...,
                    probes=300L, gdFun=mxGenerateData,
                    method=c('empirical', 'ncp'),
                    statistic=c('LRT','AIC','BIC'),
                    OK=mxOption(trueModel, "Status OK"), checkHess=FALSE,
                    silent=!interactive())
{
  garbageArguments <- list(...)
  if (length(garbageArguments) > 0) {
    stop("mxPower does not accept values for the '...' argument")
  }
  if (length(power) == 1 && is.na(power)) power <- c()
  if (length(n) == 1 && is.na(n)) n <- c()
  validateSigLevel(sig.level)
  if (length(falseModel) > 1) {
    if (length(n) > 1 || length(power) > 1) {
      stop(paste("You cannot pass more than 1 falseModel at the same time",
	"that you pass more than 1 sample size or power"))
    }
    got <- sapply(falseModel,
                  function(x) mxPower(trueModel, x, n, sig.level=sig.level,
                                      power=power, probes=probes,
                                      gdFun=gdFun, method=method, statistic=statistic, OK=OK,
                                      checkHess=checkHess))
    return(got)
  }
  if (length(n) > 1) {
    if (length(falseModel) > 1 || length(power) > 1) {
      stop(paste("You cannot pass more than 1 sample size at the same time",
	"that you pass more than 1 falseModel or power"))
    }
    got <- sapply(n,
                  function(x) mxPower(trueModel, falseModel, x, sig.level=sig.level,
                                      power=power, probes=probes,
                                      gdFun=gdFun, method=method, statistic=statistic, OK=OK,
                                      checkHess=checkHess))
    return(got)
  }
  
  origSampleSize <- totalSampleSize(trueModel)
  method <- match.arg(method)
  statistic <- match.arg(statistic)
  detail <- list(method=method, sig.level=sig.level, statistic=statistic)
  if (is.null(power)) {
    if (is.null(n)) stop("To estimate power, it is necessary to fix sample size (set n = )")
    if (method == 'ncp') {
      got <- mxPowerSearch(trueModel, falseModel, sig.level=sig.level, method=method,
                           grid=n, statistic=statistic)
      ret <- got[,'power']
      detail$n <- n
      detail$power <- ret
    } else {
      result <- data.frame(seed=as.integer(runif(probes, min = -2e9, max=2e9)),
        reject=NA, mseTrue=NA,
        statusTrue=as.statusCode(NA), statusFalse=as.statusCode(NA))
      
      trueModel <- ProcessCheckHess(trueModel, checkHess)
      falseModel <- ProcessCheckHess(falseModel, checkHess)

      prevProgressLen <- 0L

      detail$probes <- probes
      for (rx in 1:probes) {
        set.seed(result[rx,'seed'])
        info <- paste(rx, "/", probes)
        if (!silent) imxReportProgress(info, prevProgressLen)
        prevProgressLen <- nchar(info)
        simData <- try(gdFun(trueModel, returnModel=FALSE, nrowsProportion=n / origSampleSize))
        if (is(simData, "try-error")) {
          stop(paste("Cannot generate data with trueModel",
            omxQuotes(trueModel$name)), call.=FALSE)
        }
        if (is(simData, "data.frame")) {
          simData <- list(simData)
          names(simData) <- trueModel$name
        }
        
        true1  <- loadDataIntoModel(trueModel,  simData)
        true1  <- mxRun(true1,  silent=TRUE, suppressWarnings = TRUE)
        # complain about parameters at box constraints TODO
        
        topDataIndex <- match(trueModel$name, names(simData))
        names(simData)[topDataIndex] <- falseModel$name
        false1 <- loadDataIntoModel(falseModel, simData)
        false1 <- mxRun(false1, silent=TRUE, suppressWarnings = TRUE)
        if (statistic == 'LRT') {
          cmp1 <- mxCompare(true1, false1)
          pval <- cmp1[2,'p']
          if (is.na(pval)) stop("falseModel fit better than trueModel?")
          result[rx, 'reject'] <- pval < sig.level
        } else if (statistic == 'AIC') {
          rej <- summary(true1)[['informationCriteria']]['AIC:','df'] <
            summary(false1)[['informationCriteria']]['AIC:','df']
          result[rx, 'reject'] <- rej
        } else if (statistic == 'BIC') {
          rej <- summary(true1)[['informationCriteria']]['BIC:','df'] <
            summary(false1)[['informationCriteria']]['BIC:','df']
          result[rx, 'reject'] <- rej
        } else { stop(statistic) }
        result[rx, 'mseTrue'] <- sum((coef(true1) - coef(trueModel))^2)
        result[rx, 'statusTrue'] <- as.statusCode(true1$output$status$code)
        result[rx, 'statusFalse'] <- as.statusCode(false1$output$status$code)
      }
      if (!silent) imxReportProgress('', prevProgressLen)
      # return detailed results somehow? TODO
      okResult <- result[result$statusTrue %in% OK & result$statusFalse %in% OK,]
      ret <- sum(okResult$reject) / nrow(okResult)
      detail$power <- ret
      detail$n <- n
    }
  } else { # length(power) > 0
    detail$power <- power
    if (method == 'empirical') detail$probes <- probes
    if (is.null(n)) {
      # search n:power relationship
      result <- mxPowerSearch(trueModel, falseModel, probes=probes, gdFun=gdFun,
	      method=method, statistic=statistic, OK=OK, checkHess=checkHess, sig.level=sig.level)
    } else {
      # search parameter:power relationship
      result <- mxPowerSearch(trueModel, falseModel, n=n, probes=probes, gdFun=gdFun,
	      method=method, statistic=statistic, OK=OK, checkHess=checkHess, sig.level=sig.level)
	    detail$parameter <- setdiff(names(coef(trueModel)), names(coef(falseModel)))
	    detail$n <- n
    }
    model <- attr(result, "model")
    if (is.null(model)) {
      ret <- sapply(power, function(p1) {
        ind <- min(1 + findInterval(p1, result$power, all.inside = TRUE), nrow(result))
        result[ind, 'N']
      })
    } else {
      ret <- qlogis(power, -coef(model)[1]/coef(model)[2], 1/coef(model)[2])
      if (is.null(n)) ret <- round(ret * origSampleSize)
    }
    if (is.null(n)) {
	    detail$n <- ret
    } else {
	    detail$parameterDiff <- ret
    }
  }
  if (length(ret) == 1) {
	  attr(ret, 'detail') <- detail
	  class(ret) <- 'result.mxPower'
  }
  ret
}

print.result.mxPower <- function(x,...) {
	detail <- attr(x, 'detail')
	keys <- sort(names(detail))
	klen <- max(nchar(keys))
	for (kx in 1:length(keys)) {
		cat(sprintf(paste0("%", klen, "s = ",
			format(detail[[ keys[kx] ]])), keys[kx]), fill=TRUE)
	}
}
