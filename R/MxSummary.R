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

calculateConstraints <- function(model, useSubmodels) {
	constraints <- model@runstate$constraints
	retval <- c()
	if (length(constraints) > 0) {
		retval <- sapply(constraints, calculateConstraintsHelper, model)
	}
	if (useSubmodels && length(model@runstate$independents) > 0) {
		submodelConstraints <- lapply(model@runstate$independents, calculateConstraints, FALSE)
		names(submodelConstraints) <- NULL
		submodelConstraints <- unlist(submodelConstraints)
		retval <- c(retval, submodelConstraints)
	}
	return(retval)
}

calculateConstraintsHelper <- function(constraint, model) {
	if (constraint@relation == "==") {
		leftHandSide <- constraint@formula[[2]]
					# summary should not use mxEval but only examine the information in runstate
		value <- eval(substitute(mxEval(x, model, compute=TRUE),
			list(x = leftHandSide)))
		value <- as.matrix(value)
		return(nrow(value) * ncol(value))
	} else {
		return(0)
	}
}

observedStatisticsHelper <- function(model, expectation, datalist, historySet) {
	if ('numStats' %in% slotNames(expectation)) {
		if (length(expectation@numStats) > 0 && !is.na(expectation@numStats)) {
			return(list(expectation@numStats, historySet))
		}
	}
	if (length(expectation@data)==0 || is.na(expectation@data) || !.hasSlot(expectation, 'dims')) {
		return(list(0, historySet))
	}
	if (is.numeric(expectation@data)) {
		data <- datalist[[expectation@data + 1]]
	} else {
		data <- model[[expectation@data]] 
	}
	if (is(data, "MxDataDynamic")) {
		# need to revisit TODO
		return(list(0, historySet))
	}
	if (data@type == 'cov' || data@type == 'sscp') {
		if (data@name %in% historySet) {
			return (list(0, historySet))
		}
		n <- nrow(data@observed)
		dof <- n * (n + 1) / 2
		if (!single.na(data@means)) {
			dof <- dof + length(data@means)
		}
		historySet <- append(data, historySet)
	} else if (data@type == 'cor') {
		if (data@name %in% historySet) {
			return (list(0, historySet))
		}
		n <- nrow(data@observed)
		dof <- n * (n - 1) / 2
		if (!single.na(data@means)) {
			dof <- dof + length(data@means) 
		}
		historySet <- append(data, historySet)
	} else if (is(expectation, "MxExpectationBA81")) {  # refactor TODO
		if (is.na(expectation@weightColumn)) return(list(NA, historySet))
		dof <- nrow(data@observed) - 1
		historySet <- append(data, historySet)
	} else if (data@type == 'acov') {
		if (data@name %in% historySet) {
			return (list(0, historySet))
		}
		numThresh <- sum(!is.na(data@thresholds))
		numMeans <- sum(!is.na(data@means))
		n <- nrow(data@observed)
		# Include diagonal of observed when no thresholds
		dof <- numThresh + numMeans + ifelse(numThresh > 0, n*(n-1)/2, n*(n+1)/2)
		historySet <- append(data, historySet)
	} else {
		dof <- 0
		observed <- data@observed
		for (i in 1:ncol(observed)) {
			colname <- colnames(observed)[[i]]
			fullname <- paste(data@name, colname, sep='.')
			if ((colname %in% expectation@dims) && !(fullname %in% historySet)) {
				dof <- dof + sum(!is.na(observed[,i]))
				historySet <- append(fullname, historySet)
			}
		}
	}
	return(list(dof, historySet))
}

observedStatistics <- function(model, useSubmodels, constraintOffset) {
	datalist <- model@runstate$datalist
	expectations <- model@runstate$expectations
	retval <- constraintOffset
	if (length(expectations) > 0) {
		historySet <- character()
		for(i in 1:length(expectations)) {
			result <- observedStatisticsHelper(model, expectations[[i]], datalist, historySet)
			retval <- retval + result[[1]]
			historySet <- result[[2]]
		}
	}
	if (useSubmodels && length(model@runstate$independents) > 0) {
		submodelStatistics <- sapply(model@runstate$independents, observedStatistics, FALSE, 0)
		retval <- retval + sum(submodelStatistics)
	}
	return(retval)
}

fitfunctionNumberObservations <- function(fitfunction) {
	if (("numObs" %in% slotNames(fitfunction)) && !single.na(fitfunction@numObs)) {
		return(fitfunction@numObs)
	} else {
		return(0)
	}
}

numberObservations <- function(datalist, fitfunctions) {
	datalist <- Filter(function(x) !is(x,"MxDataDynamic"), datalist)
	dataObservations <- sapply(datalist, slot, name = "numObs")
	fitfunctionObservations <- sapply(fitfunctions, fitfunctionNumberObservations)
	return(sum(as.numeric(dataObservations), as.numeric(fitfunctionObservations)))
}

computeFValue <- function(datalist, likelihood, chi) {
	if(length(datalist) == 0) return(NA)
	datalist <- Filter(function(x) !is(x,"MxDataDynamic"), datalist)
	if(all(sapply(datalist, function(x) 
		{x@type == 'raw'}))) return(likelihood)
	if(all(sapply(datalist, function(x) 
		{x@type == 'cov'}))) return(chi)
	if(all(sapply(datalist, function(x) 
		{x@type == 'cor'}))) return(chi)
	if(all(sapply(datalist, function(x) 
		{x@type == 'acov'}))) return(chi)
	return(NA)
}

computeFitStatistics <- function(likelihood, DoF, chi, chiDoF, numObs,
				 independence, indDoF, saturated=0, satDoF=0) {
	if (!is.na(independence) && !is.na(likelihood) && independence < likelihood) {
		warning(paste("Your model may be mis-specified (and fit worse than an independence model),",
			      "or you may be using the wrong independence model, see ?mxRefModels"))
	}
	CFI <- (independence - indDoF - likelihood + DoF)/(independence - indDoF - saturated + satDoF)
	TLI <- 1
	rmseaSquared <- 0
	RMSEA <- 0
	RMSEACI <- c(lower=NA, upper=NA)
	RMSEANull <- 0.05
	RMSEAClose <- NA
	if (!is.na(chiDoF) && chiDoF > 0) {
		TLI <- ((independence-saturated)/(indDoF-satDoF) - (chi)/(DoF-satDoF))/((independence-saturated)/(indDoF-satDoF) - 1)
					# Here we use N in the denominator as given in the original
					# RMSEA paper. The difference between N and N-1 is negligible
					# for sample sizes over 30. RMSEA should not be taken seriously
					# such small samples anyway.
		rmseaSquared <- (chi / (chiDoF) - 1) / numObs
		if (length(rmseaSquared) == 0 || is.na(rmseaSquared) || 
		    is.nan(rmseaSquared)) { 
					# || (rmseaSquared < 0)) { # changed so 'rmseaSquared < 0' yields zero with comment
			RMSEA <- NA
		} else {
			RMSEA <- ifelse(rmseaSquared < 0, 0.0, sqrt(rmseaSquared))
			ci <- try(rmseaConfidenceIntervalHelper(chi, chiDoF, numObs, .025, .975))
			if (!inherits(ci, "try-error")) RMSEACI <- ci
			RMSEAClose <- rmseaPCloseHelper(chi, chiDoF, numObs, null=RMSEANull)
		}
	}
	list(CFI=CFI, TLI=TLI, RMSEA=RMSEA, RMSEASquared=rmseaSquared, RMSEACI=RMSEACI, RMSEAClose=RMSEAClose, RMSEANull=RMSEANull)
}

catFitStatistics <- function(x) {
	cat("CFI:", x$CFI, '\n')
	cat("TLI:", x$TLI, '  (also known as NNFI)', '\n')
	if (length(x$RMSEASquared) == 1 && !is.na(x$RMSEASquared) && x$RMSEASquared < 0.0) {
	  cat("RMSEA: ", x$RMSEA, " *(Non-centrality parameter is negative)", "  [95% CI (", x$RMSEACI[1], ", ", x$RMSEACI[2], ")]", '\n', sep="")
	} else {
		cat("RMSEA:  ", x$RMSEA, "  [95% CI (", x$RMSEACI[1], ", ", x$RMSEACI[2], ")]", '\n', sep="")
	}
	cat("Prob(RMSEA <= ", x$RMSEANull, "): ", x$RMSEAClose, '\n', sep='')
}

fitStatistics <- function(model, useSubmodels, retval) {
	datalist <- model@runstate$datalist
	likelihood <- retval[['Minus2LogLikelihood']]
	saturated <- retval[['SaturatedLikelihood']]
	independence <- retval[['IndependenceLikelihood']]
	if(is.null(model@output$chi)){
		chi <- likelihood - saturated
	} else {chi <- model@output$chi}
	DoF <- retval$degreesOfFreedom
	satDoF <- retval$saturatedDoF
	indDoF <- retval$independenceDoF
	nParam <- dim(retval$parameters)[1]
	Fvalue <- computeFValue(datalist, likelihood, chi)
	if(is.null(model@output$chiDoF)){
		chiDoF <- DoF - satDoF # DoF = obsStat-model.ep; satDoF = obsStat-sat.ep; So sat.ep-model.ep == DoF-satDoF
	} else {chiDoF <- model@output$chiDoF}
	retval[['ChiDoF']] <- chiDoF
	retval[['Chi']] <- chi
	retval[['p']] <- suppressWarnings(pchisq(chi, chiDoF, lower.tail = FALSE))
	retval[['AIC.Mx']] <- Fvalue - 2 * DoF
	retval[['BIC.Mx']] <- (Fvalue - DoF * log(retval[['numObs']])) 
	AIC.p <- Fvalue + 2 * nParam
	BIC.p <- (Fvalue + nParam * log(retval[['numObs']])) 
	sBIC <- (Fvalue + nParam * log((retval[['numObs']]+2)/24))
	retval[['satDoF']] <- satDoF
	retval[['indDoF']] <- indDoF
	IC <- matrix(NA, nrow=2, ncol=3, dimnames=list(c("AIC:", "BIC:"), c('df', 'par', 'sample')))
	IC[,'df'] <- c(retval$AIC.Mx, retval$BIC.Mx)
	IC[,'par'] <- c(AIC.p, BIC.p)
	IC['BIC:','sample'] <- sBIC
	retval[['informationCriteria']] <- IC
	
	retval$fitUnits <- model@output$fitUnits
	
	fi <- computeFitStatistics(likelihood, DoF, chi, chiDoF,
		retval[['numObs']], independence, indDoF, saturated, satDoF)
	for (k in names(fi)) retval[[k]] <- fi[[k]]
	return(retval)
}

omxRMSEA <- function(model, lower=.025, upper=.975, null=.05, ...){
	smod <- summary(model, ...)
	x2 <- smod$Chi
	df <- smod$ChiDoF
	N <- smod$numObs
	rmsea <- smod$RMSEA
	ci <- rmseaConfidenceIntervalHelper(chi.squared=x2, df=df, N=N, lower=lower, upper=upper)
	pn <- rmseaPCloseHelper(x2, df, N, null)
	return(c(ci[1], est.rmsea=rmsea, ci[2], null=null, `Prob(x <= null)`=pn))
}

rmseaConfidenceIntervalHelper <- function(chi.squared, df, N, lower, upper){
	# Lower confidence interval
	if( pchisq(chi.squared, df=df, ncp=0) >= upper){ #sic
		lower.lam <- uniroot(f=pChiSqFun, interval=c(1e-10, 1e4), val=chi.squared, degf=df, goal=upper)$root
		# solve pchisq(ch, df=df, ncp=x) == upper for x
	} else{
		lower.lam <- 0
	}
	# Upper confidence interval
	if( pchisq(chi.squared, df=df, ncp=0) >= lower){ #sic
		upper.lam <- uniroot(f=pChiSqFun, interval=c(1e-10, 1e4), val=chi.squared, degf=df, goal=lower)$root
		# solve pchisq(ch, df=df, ncp=x) == lower for x
	} else{
		upper.lam <- 0
	}
	lower.rmsea <- sqrt(lower.lam/(N*df))
	upper.rmsea <- sqrt(upper.lam/(N*df))
	return(c(lower=lower.rmsea, upper=upper.rmsea))
}

pChiSqFun <- function(x, val, degf, goal){
	goal - pchisq(val, degf, ncp=x)
}

rmseaPCloseHelper <- function(x2, df, N, null){
	1-pchisq(x2, df=df, ncp=N*df*(null)^2)
}


parameterList <- function(model, useSubmodels) {
	if (useSubmodels && length(model@runstate$independents) > 0) {
		ptable <- parameterListHelper(model, TRUE)
		submodelParameters <- lapply(model@runstate$independents, parameterListHelper, TRUE)
		ptable <- Reduce(rbind, submodelParameters, ptable)
	} else {
		ptable <- parameterListHelper(model, FALSE)
	}
	return(ptable)
}

parameterListHelper <- function(model, withModelName) {
	ptable <- data.frame()
	if(length(model@output) == 0) { return(ptable) }
	estimates <- model@output$estimate
    if (!is.null(model@output$standardErrors) && 
 		length(model@output$standardErrors) == length(estimates)) { 
 		errorEstimates <- model@output$standardErrors 
	} else { 
		errorEstimates <- rep.int(as.numeric(NA), length(estimates)) 
	}
	matrices <- model@runstate$matrices
	parameters <- model@runstate$parameters
	if (length(estimates) > 0) {
		matrixNames <- names(matrices)
		for(i in 1:length(estimates)) {
			mLocation <- parameters[[i]][[5]][[1]] + 1
			mRow <- parameters[[i]][[5]][[2]] + 1
			mCol <- parameters[[i]][[5]][[3]] + 1
			lbound <- parameters[[i]][[1]]
			ubound <- parameters[[i]][[2]]
			aMatrix <- matrices[[mLocation]][[1]]
			if (getOption('mxShowDimnames')) {
				if (!is.null(rownames(aMatrix))) {
					mRow <- rownames(aMatrix)[[mRow]]
				}
				if (!is.null(colnames(aMatrix))) {
					mCol <- colnames(aMatrix)[[mCol]]
				}
			}
			if (withModelName) {
				ptable[i, 'model'] <- model@name
			}
			if (!is.null(names(estimates))) {
				ptable[i, 'name'] <- names(estimates)[[i]]
			} else {
				ptable[i, 'name'] <- paste("p", i, sep="")
			}
			ptable[i, 'matrix'] <- simplifyName(matrixNames[[mLocation]], model@name)
			ptable[i, 'row'] <- mRow
			ptable[i, 'col'] <- mCol
			ptable[i, 'Estimate'] <- estimates[[i]]
			ptable[i, 'Std.Error'] <- errorEstimates[[i]]
			ptable[i, 'lbound'] <- lbound
			ptable[i, 'ubound'] <- ubound
		}
	}
	return(ptable)
}

computeOptimizationStatistics <- function(model, numStats, useSubmodels, saturatedDoF, independenceDoF, retval) {
	# get estimated parameters
	estimates <- model@output$estimate
	# should saturated/independence models include means?
	if(length(model@runstate$datalist)==1){
		type <- model@runstate$datalist[[1]]@type
		means <- model@runstate$datalist[[1]]@means
		# if there's raw data, then use means in saturated/independence models
		if(type=="raw"){
			useMeans <- TRUE
		} else {
		# if there's not raw data, only use means if they're present
			if((dim(means)[2]==1)&is.na(means[1,1])){
				useMeans <- FALSE
			} else{
				useMeans <- TRUE	
			}
		}
		# number of variables
		if(model@runstate$datalist[[1]]@type != 'raw'){
			nvar <- dim(model@runstate$datalist[[1]]@observed)[2]
		} else if( length(model@runstate$expectations) == 1 ) {
			nvar <- length(model@runstate$expectations[[1]]@dims)
		} else {
			nvar <- 0
		}
	# if there are multiple or zero datalists, then do nothing
	} else {
		useMeans <- NA	
		nvar <- 0
	}
		# how many thresholds does each variable have (needed for saturated and independence DoF calculation)
	# grab the expectation
	obj <- model@runstate$expectation
	# grab the thresholdLevels object and expected means; punt if there is more than one expectation
	if (length(obj)==1){
		if ("thresholdLevels" %in% slotNames(obj[[1]])){
			thresholdLevels <- obj[[1]]@thresholdLevels
			if (length(thresholdLevels)==0){thresholdLevels <- rep(NA, nvar)}
		} else {
			thresholdLevels <- rep(NA, nvar)
		}
	} else {
		thresholdLevels <- NULL	
	}
	# number of continuous variables, provided there is just one expectation
	if (!is.null(thresholdLevels)){
		continuous <- sum(is.na(thresholdLevels))
	} else{
		continuous <- NA
	}
	# number of thresholds in the model
	if (!is.null(thresholdLevels)){
		thresh <- sum(thresholdLevels, na.rm=TRUE)
	} else{
		thresh <- NA
	}
	# constraints, parameters, model degrees of freedom
	retval[['constraints']] <- calculateConstraints(model, useSubmodels)
	retval[['estimatedParameters']] <- nrow(retval$parameters)
  if(any(sapply(obj,function(x){"MxExpectationGREML" %in% class(x)}))){
    retval[['estimatedParameters']] <- retval[['estimatedParameters']] + 
      sum(sapply(obj,imxExtractSlot,name="numFixEff"))
  }
	if (is.null(numStats)) {
		retval[['observedStatistics']] <- observedStatistics(model, useSubmodels, sum(retval$constraints))
	} else {
		retval[['observedStatistics']] <- numStats
	}
	retval[['degreesOfFreedom']] <- retval$observedStatistics - retval$estimatedParameters
	# calculate or populate saturated degrees of freedom
	if(is.null(saturatedDoF)) {
		retval[['saturatedDoF']] <- retval$observedStatistics - (nvar * (nvar-1) / 2 + continuous*(1+useMeans) + thresh)
	} else {
		retval[['saturatedDoF']] <- saturatedDoF
	}
	# calculate or populate independence degrees of freedom
	if(is.null(independenceDoF)) {
		# indDoF = 1 df per continuous variable variance + 1 df per continuous mean + 1 df per threshold
		retval[['independenceDoF']] <- retval$observedStatistics - (continuous*(1+useMeans) + thresh)
	} else {
		retval[['independenceDoF']] <- independenceDoF
	}
	# set NULLs to NAs
	if (is.null(retval$saturatedDoF)) {
		retval$SaturatedDoF <- NA
	}
	if (is.null(retval$independenceDoF)) {
		retval$IndependenceDoF <- NA
	}
	retval[['saturatedParameters']] <- retval[['observedStatistics']] - retval[['saturatedDoF']]
	retval[['independenceParameters']] <- retval[['observedStatistics']] - retval[['independenceDoF']]
	# calculate fit statistics
	retval <- fitStatistics(model, useSubmodels, retval)
	return(retval)
}

print.summary.mxmodel <- function(x,...) {
	cat("Summary of", x$modelName, '\n', '\n')
	if(x$verbose==TRUE){
		if (length(x$compute)) {
			cat("compute plan:\n")
			print(x$compute)
			cat("\n")
		}
		if (length(x$dataSummary) > 0) {
			cat("data:\n")
			print(x$dataSummary)
		}
	}
	if (!is.null(x$npsolMessage)) {
		cat(x$npsolMessage,'\n','\n')
	}
	params <- x$parameters
	if (!is.null(params) && nrow(params)) {
		cat("free parameters:\n")
		params$lbound <- mapply(highlightProblem, params$lbound, params$lboundMet)
		params$ubound <- mapply(highlightProblem, params$ubound, params$uboundMet)
		params$lbound[is.na(params$lbound)] <- ""
		params$ubound[is.na(params$ubound)] <- ""
		params$lboundMet <- NULL
		params$uboundMet <- NULL
		if (!x$verbose) {
			if (all(is.na(params[['Std.Error']]))) {
				params[['Std.Error']] <- NULL
			}
			if (all(params$lbound == "" & params$ubound == "")) {
				params[['lbound']] <- NULL
				params[['ubound']] <- NULL
			}
		}
		if (!is.null(params[['Std.Error']]) && !is.null(x[['seSuspect']])) {
			seCol <- match('Std.Error', colnames(params))
			before <- params[1:seCol]
			stars <- mapply(highlightProblem, rep('',length(x[['seSuspect']])), x[['seSuspect']])
			if (length(params) > seCol) {
				params <- cbind(before, 'A'=stars, params[(seCol+1):length(params)])
			} else {
				params <- cbind(before, 'A'=stars)
			}
		}
		if (!is.null(x$bootstrapQuantile) && nrow(x$bootstrapQuantile) == nrow(params)) {
			bq <- x$bootstrapQuantile
			params <- cbind(params, bq)
		}
		cmap <- 1:ncol(params)
		isBound <- colnames(params) %in% paste0(c('l','u'),'bound')
		params <- params[,c(cmap[!isBound], cmap[isBound]),drop=FALSE]
		print(params)
		cat('\n')
	}
	if (!is.null(x$CI) && length(x$CI) > 0) {
		cat("confidence intervals:\n")
		print(x$CI)
		if (any(is.na(x$CI[,c('lbound','ubound')])) && !(x$verbose)) {
			cat("  To investigate missing CIs, run summary() again, with verbose=T, to see CI details.", '\n')
		}
		cat('\n')
	}
	if(x$verbose && length(x$CIdetail)){
		cat("CI details:\n")
		print(x$CIdetail)
		cat("\n")
	}
	if(length(x$GREMLfixeff)>0 && any(sapply(x$GREMLfixeff,length)>0)){
		cat("regression coefficients:\n")
		print(x$GREMLfixeff)
		cat("\n")
	}
	cat('Model Statistics:', '\n')
	EP <- matrix(
		c(x$estimatedParameters, x$saturatedParameters, x$independenceParameters,
		x$degreesOfFreedom, x$saturatedDoF, x$independenceDoF,
		x$Minus2LogLikelihood, x$SaturatedLikelihood, x$IndependenceLikelihood),
		nrow=3, ncol=3,
		dimnames=list(
			c('       Model:', '   Saturated:', 'Independence:'),
			c(' |  Parameters', ' |  Degrees of Freedom', paste0(' |  Fit (', x$fitUnits, ' units)'))
		)
	)
	print(EP)
	cat('Number of observations/statistics: ', x$numObs, "/", x$observedStatistics, '\n\n', sep="")
	constraints <- x$constraints
	if(length(constraints) > 0) {
		for(i in 1:length(constraints)) {
			name <- names(constraints)[[i]]
			if (constraints[[i]] == 1) plural <- ''
			else plural <- 's'
			cat("Constraint", omxQuotes(simplifyName(name, x$modelName)), "contributes",
				constraints[[i]], paste("observed statistic", plural, '.', sep=''), "\n")
		}
	}
	if (!is.null(x$infoDefinite) && !is.na(x$infoDefinite)) {
		if (!x$infoDefinite) {
			cat("\n** Information matrix is not positive definite (not at a candidate optimum).\n  Be suspicious of these results. At minimum, do not trust the standard errors.\n\n")
		} else if(x$verbose==TRUE) {
			cat("condition number of the information matrix: ", x$conditionNumber, "\n")
		}
	}
	if (x$verbose==TRUE && !is.null(x$maxAbsGradient)) {
		cat("maximum absolute gradient: ", x$maxAbsGradient, " (",names(x$maxAbsGradient),")\n")
	}
	#
	# Chi-square goodness of fit test
	if(x$verbose==TRUE || !is.na(x$Chi)){
		chival <- x$Chi
		if(is.na(x$SaturatedLikelihood)){
			chidof <- NA
		} else {
			chidof <- x$ChiDoF
		}
		chipee <- x$p
		cat("chi-square:  ", "\U03C7\U00B2 ( df=", chidof, " ) = ", chival, ",  p = ", chipee, '\n', sep="")
	}
	#
	cat("Information Criteria: \n")
	IC <- x$informationCriteria
	colnames(IC) <- c(" |  df Penalty", " |  Parameters Penalty", " |  Sample-Size Adjusted")
	print(IC)
	#
	# Absolute fit indices
	if(x$verbose==TRUE || any(!is.na(c(x$CFI, x$TLI, x$RMSEA)))){
		catFitStatistics(x)
	}
	if(any(is.na(c(x$CFI, x$TLI, x$RMSEA)))){
		cat("To get additional fit indices, see help(mxRefModels)\n")
	}
	#
	# Timing information
	cat("timestamp:", format(x$timestamp), '\n')
	if(x$verbose==TRUE){
		cat("frontend time:", format(x$frontendTime), '\n')
		cat("backend time:", format(x$backendTime), '\n')
		cat("independent submodels time:", format(x$independentTime), '\n')
		cat("cpu time:", format(x$cpuTime), '\n')
	}
	cat("Wall clock time (HH:MM:SS.hh):", format(.POSIXct(x$wallTime, tz="GMT"), "%H:%M:%OS2"), "\n")
	if (x$verbose==FALSE && !is.null(x$optimizerEngine)) {
		cat("optimizer: ", x$optimizerEngine, '\n')
	}
	cat("OpenMx version number:", format(x$mxVersion), '\n')
	cat("Need help?  See help(mxSummary)", '\n')
	cat('\n')
	if (!x$wasRun) {
		message("WARNING: This model has not been run yet. Tip: Use\n  model = mxRun(model)\nto estimate a model.")
	} else if (x$stale) {
		message("WARNING: This model was modified since it was run. Summary information may be out-of-date.")
	}
}

setLikelihoods <- function(model, saturatedLikelihood, independenceLikelihood, retval) {
	# populate saturated -2 log likelihood
	if(is.null(saturatedLikelihood)) {
		retval$SaturatedLikelihood <- attr(model@fitfunction$result, "SaturatedLikelihood")
	} else {
		retval$SaturatedLikelihood <- saturatedLikelihood
	}
	# populate independence -2 log likelihood	
	if(is.null(independenceLikelihood)) {
		retval$IndependenceLikelihood <- attr(model@fitfunction$result, "IndependenceLikelihood")
	} else {
		retval$IndependenceLikelihood <- independenceLikelihood
	}
	# populate model -2 log likelihood
	retval$Minus2LogLikelihood <- model@output$Minus2LogLikelihood
	# set NULLs to NAs
	if (is.null(retval$SaturatedLikelihood)) {
		retval$SaturatedLikelihood <- NA
	}
	if (is.null(retval$Minus2LogLikelihood)) {
		retval$Minus2LogLikelihood <- NA
	}
	if (is.null(retval$IndependenceLikelihood)) {
		retval$IndependenceLikelihood <- NA
	}
	return(retval)
}

setNumberObservations <- function(numObs, datalist, fitfunctions, retval) {
	if(is.null(numObs)) {
		retval$numObs <- numberObservations(datalist, fitfunctions)
	} else {
		retval$numObs <- numObs
	}
	return(retval)
}

generateDataSummary <- function(model, useSubmodels) {
	datalist <- model@runstate$datalist
	retval <- lapply(datalist, summarize)
	if (useSubmodels && length(model@runstate$independents) > 0) {
		submodelSummary <- lapply(model@runstate$independents, generateDataSummary, FALSE)
		names(submodelSummary) <- NULL
		submodelSummary <- unlist(submodelSummary, recursive = FALSE)
		retval <- c(retval, submodelSummary)
	}
	return(retval)
}

##' imxEvalByName
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' This function should not be used in MxSummary. All summary
##' information should be extracted from runstate.
##'
##' @param name name
##' @param model model
##' @param compute compute
##' @param show show
imxEvalByName <- function(name, model, compute=FALSE, show=FALSE) {
   if ((length(name) != 1) || typeof(name) != "character") {
      stop("'name' argument must be a character argument")
   }
   if (!is(model, "MxModel")) {
      stop("'model' argument must be a MxModel object")
   }
   if (hasSquareBrackets(name)) {
      components <- splitSubstitution(name)
      eval(substitute(mxEval(x[y,z], model, compute, show),
         list(x = as.name(components[[1]]), 
            y = parse(text=components[[2]])[[1]],
            z = parse(text=components[[3]])[[1]])))
   } else {
      eval(substitute(mxEval(x, model, compute, show),
         list(x = as.name(name))))
   }
}

boundsMet <- function(model, retval){
	params <- retval$parameters
	lbound <- params$lbound
	ubound <- params$ubound
	estimate <- params$Estimate
	threshold <- model@options[["Feasibility tolerance"]]
	if (is.null(threshold)){
		threshold <- getOption("mxOptions")[["Feasibility tolerance"]]
	}
	threshold <- as.numeric(threshold)
	lboundMet <- mapply(compareBounds, estimate, lbound, MoreArgs=list(threshold))
	uboundMet <- mapply(compareBounds, estimate, ubound, MoreArgs=list(threshold))
	params$lboundMet <- lboundMet
	params$uboundMet <- uboundMet
	retval$parameters <- params
	return(retval)
}

compareBounds <- function(estimate, bound, threshold){
	if (is.na(bound)){
		return(FALSE)
	}
	absDelta <- abs(estimate - bound)
	return (absDelta < threshold)
}

highlightProblem <- function(bound, boundMet){
	if (boundMet){
		if (is.numeric(bound)) bound <- round(bound,4)
		return(paste(bound, "!", sep=""))
	}
	else {
		return(bound)
	}
}

parseLikelihoodArg <- function(input, arg) {
	input <- input[[arg]]
	if (is.null(input)) {
		return(input)
	} else if (is.numeric(input)) {
		return(input)
	} else if (is(input, "MxModel")) {
		name <- input@name
		if (is.null(input@fitfunction)) {
			stop(paste(omxQuotes(name), "model passed",
				"to summary function does not",
				"have top-level fitfunction in",
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (length(input@fitfunction@result) != 1) {
			stop(paste(omxQuotes(name), "model passed to summary",
				"function does not have a 1x1 matrix",
				"result in fitfunction in",
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		return(input@fitfunction@result[1,1])
	} else if(is.list(input) && length(input)==2) {
		stop(paste("List of length two (illegal argument) passed to", omxQuotes(arg),
			   "argument of summary function. You probably meant to use",
			   "the refModels argument instead."), call. = FALSE)
	} else {
		stop(paste("Illegal argument passed to", omxQuotes(arg),
			"argument of summary function in",
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
}

parseDfArg <- function(input, arg) {
	input <- input[[arg]]
	if (is.null(input)) {
		return(input)
	} else if (is.numeric(input)) {
		return(input)
	} else if (is(input, "MxModel")) {
		name <- input@name
		if (is.null(input@fitfunction)) {
			stop(paste(omxQuotes(name), "model passed",
				"to summary function does not",
				"have top-level fitfunction in",
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		if (length(input@fitfunction@result) != 1) {
			stop(paste(omxQuotes(name), "model passed to summary",
				"function does not have a 1x1 matrix",
				"result in fitfunction in",
				deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
		}
		return(summary(input)$degreesOfFreedom)
	} else {
		stop(paste("Illegal argument passed to", omxQuotes(arg),
			"argument of summary function in",
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
}

refToLikelihood <- function(model) {
	if (is(model, "MxModel")) {
		if (!model@.wasRun) stop("Reference model must be run to obtain fit indices")
		model$output$Minus2LogLikelihood
	} else if (is.list(model)) {
		model[[1]]
	} else {
		stop(paste("Illegal argument passed to refModels",
			"argument of summary function in",
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
}

refToDof <- function(model) {
	if (is(model, "MxModel")) {
		if (!model@.wasRun) stop("Reference model must be run to obtain fit indices")
		return(summary(model)$degreesOfFreedom)
	} else if (is.list(model)) {
		model[[2]]
	} else {
		stop(paste("Illegal argument passed to refModels",
			"argument of summary function in",
			deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
	}
}

summary.MxModel <- function(object, ..., verbose=FALSE) {
	model <- object
	dotArguments <- list(...)
	if (!is.null(dotArguments[["refModels"]])) {
		refModels <- dotArguments[["refModels"]]
		satModel <- refModels[['Saturated']]
		indModel <- refModels[['Independence']]
		saturatedLikelihood <- refToLikelihood(satModel)
		saturatedDoF <- refToDof(satModel)
		independenceLikelihood <- refToLikelihood(indModel)
		independenceDoF <- refToDof(indModel)
	} else {
		saturatedLikelihood <- parseLikelihoodArg(dotArguments, "SaturatedLikelihood")
		saturatedDoF <- parseDfArg(dotArguments, "SaturatedDoF")
		independenceLikelihood <- parseLikelihoodArg(dotArguments, "IndependenceLikelihood")
		independenceDoF <- parseDfArg(dotArguments, "IndependenceDoF")
	}
	numObs <- dotArguments$numObs
	numStats <- dotArguments$numStats
	useSubmodels <- dotArguments$indep
	if (is.null(useSubmodels)) { useSubmodels <- TRUE }
	retval <- list(wasRun=model@.wasRun, stale=model@.modifiedSinceRun)
	retval$parameters <- parameterList(model, useSubmodels)
	if (!is.null(model@compute$steps[['ND']]) && model@compute$steps[['ND']]$checkGradient &&
	    !is.null(model@compute$steps[['ND']]$output$gradient)) {
		retval$seSuspect <- !model@compute$steps[['ND']]$output$gradient[,'symmetric']
	}
	if (is(model@compute, "MxComputeBootstrap")) {
		retval$bootstrapQuantile <- model@compute$output$result
	}
	retval$GREMLfixeff <- GREMLFixEffList(model)
	retval$infoDefinite <- model@output$infoDefinite
	retval$conditionNumber <- model@output$conditionNumber
	if (length(model@output$gradient)) {
		agrad <- abs(model@output$gradient)
		retval$maxAbsGradient <- agrad[ order(-agrad)[1] ]
	}
	retval <- boundsMet(model, retval)
	retval <- setLikelihoods(model, saturatedLikelihood, independenceLikelihood, retval)
	retval <- setNumberObservations(numObs, model@runstate$datalist, model@runstate$fitfunctions, retval)
	retval <- computeOptimizationStatistics(model, numStats, useSubmodels, saturatedDoF, independenceDoF, retval)
	retval$dataSummary <- generateDataSummary(model, useSubmodels)
	retval$CI <- as.data.frame(model@output$confidenceIntervals)
	if (length(retval$CI) && nrow(retval$CI)) {
		retval$CI <- cbind(retval$CI, note=apply(retval$CI, 1, function(ci) {
					# This should probably take into account whether both bounds
					# were requested and consider the optimizer codes also. TODO
			if (any(is.na(ci)) || ci[1] == ci[3] || ci[1] >= ci[2] || ci[2] >= ci[3]) {
				"!!!"
			} else {
				""
			}
		}))
	}
	retval$CIcodes <- model@output$confidenceIntervalCodes
	if (!is.null(model@output$status$code)) {
		message <- optimizerMessages[[as.character(model@output$status$code)]]
		retval[['npsolMessage']] <- message
	}
	if( .hasSlot(model,"compute") && length(model$compute$steps$CI) ){
		retval$CIdetail <- model$compute$steps$CI$output$detail
	}
	retval$timestamp <- model@output$timestamp
	retval$frontendTime <- model@output$frontendTime
	retval$backendTime <- model@output$backendTime
	retval$independentTime <- model@output$independentTime
	retval$wallTime <- model@output$wallTime
	retval$cpuTime <- model@output$cpuTime
	retval$mxVersion <- model@output$mxVersion
	retval$modelName <- model@name
	plan <- model@runstate$compute
	if (is(plan, "MxComputeSequence")) {
		gd <- plan$steps[['GD']]
		if (is(gd, "MxComputeGradientDescent")) {
			retval$optimizerEngine <- gd$engine
		}
	}
	retval$verbose <- verbose
	class(retval) <- "summary.mxmodel"
	return(retval)
}

logLik.MxModel <- function(object, ...) {
	model <- object
	if (model@.wasRun && model@.modifiedSinceRun) {
		msg <- paste("MxModel", omxQuotes(model@name), "was modified",
			     "since it was run. The log likelihood may be out-of-date.")
		warning(msg)
	}
	ll <- NA
	if (!is.null(model@output) & !is.null(model@output$Minus2LogLikelihood)) {
		ll <- -0.5*model@output$Minus2LogLikelihood
	}

	nobs <- numberObservations(model@runstate$datalist, model@runstate$fitfunctions)
	if (!is.na(nobs)) {
		attr(ll,"nobs") <- nobs
	}

	if (!is.null(model@output))
		attr(ll,"df")<- length(model@output$estimate)
	else
		attr(ll,"df") <- NA
	class(ll) <- "logLik"
	return(ll)
}


.standardizeParams <- function(x=NULL, model, Apos, Spos, give.matrices=FALSE){
  if(is.null(x)){x <- omxGetParameters(model)}
  param <- omxGetParameters(model)
  paramNames <- names(param)
  model <- omxSetParameters(model, values=x, labels=paramNames, free=TRUE)
  model_A <- model[[model$expectation$A]] #<--the A matrix might not be named "A".
  A <- list( model_A$values, model_A$result )
  A <- A[[which.max(c( length(A[[1]]), length(A[[2]]) ))]]
  model_S <- model[[model$expectation$S]] #<--Likewise for S
  S <- list( model_S$values, model_S$result )
  S <- S[[which.max(c( length(S[[1]]), length(S[[2]]) ))]]
  I <- diag(1, nrow(A))
  ImAInv <- solve(I-A)
  SD <- sqrt(diag(ImAInv %*% S %*% t(ImAInv)))
  SD <- diag(SD,nrow=length(SD)) #<--Needed if line immediately above is a scalar.
  InvSD <- 1/diag(SD)
  InvSD <- diag(InvSD,nrow=length(InvSD))
  Az <- InvSD %*% A %*% SD
  Sz <- InvSD %*% S %*% InvSD
  sparam <- c(Az[!is.na(Apos)],Sz[!is.na(Spos)])
  names(sparam) <- c(Apos[!is.na(Apos)],Spos[!is.na(Spos)])
  if(!give.matrices){return(sparam)}
  else{return(list(sparam=sparam,Az=Az,Sz=Sz))}
}
.mxStandardizeRAMhelper <- function(model,SE=FALSE,ParamsCov,inde.subs.flag=FALSE){
  #Recur the function for the appropriate submodels, if any:
  if(length(model@submodels)>0){
    return(lapply(
      model@submodels[which(
        sapply(model@submodels,function(x){class(x$expectation)})=="MxExpectationRAM" | 
          sapply(model@submodels,function(x){length(x@submodels)>0})  
      )],
      .mxStandardizeRAMhelper,SE=SE,ParamsCov=ParamsCov,inde.subs.flag=inde.subs.flag))
  }
  #Get A and S:
  model_A <- model[[model$expectation$A]] #<--Necessary because the A matrix might not be named "A".
  A <- list( model_A$values, model_A$result )
  A <- A[[which.max(c( length(A[[1]]), length(A[[2]]) ))]]
  model_S <- model[[model$expectation$S]] #<--Likewise for S
  S <- list( model_S$values, model_S$result )
  S <- S[[which.max(c( length(S[[1]]), length(S[[2]]) ))]]
  #Find positions of nonzero paths:
  Apos <- matrix(NA,nrow=nrow(A),ncol=ncol(A),dimnames=dimnames(A))
  Spos <- matrix(NA,nrow=nrow(S),ncol=ncol(S),dimnames=dimnames(S))
  A_need_pos <- which(A!=0,arr.ind=T)
  S_need_pos <- which(S!=0,arr.ind=T)
  S_need_pos <- subset(S_need_pos, S_need_pos[,1]>=S_need_pos[,2]) #<--Lower tri only
  numelem <- nrow(A_need_pos)+nrow(S_need_pos)
  #Create output object:
  out <- data.frame(name=vector(mode="character",length=numelem),label=vector(mode="character",length=numelem),
                    matrix=vector(mode="character",length=numelem),
                     row=vector(mode="character",length=numelem),col=vector(mode="character",length=numelem),
                    Raw.Value=vector(mode="numeric",length=numelem),
                    Raw.SE=vector(mode="numeric",length=numelem),
                    Std.Value=vector(mode="numeric",length=numelem),
                    Std.SE=vector(mode="numeric",length=numelem),stringsAsFactors=FALSE)
  out$label <- NA
  j <- 1
  #Create position strings where needed and begin to populate output:
  if(nrow(A_need_pos)>0){
    for(i in 1:nrow(A_need_pos)){
      Apos[A_need_pos[i,1],A_need_pos[i,2]] <- out$name[j] <- paste(
        model@name,".A[",A_need_pos[i,1],",",A_need_pos[i,2],"]",sep="")
      if(!all.na(model_A$labels)){
        out$label[j] <- model_A$labels[A_need_pos[i,1],A_need_pos[i,2]]
      }
      out$matrix[j] <- "A"
      out$Raw.Value[j] <- A[A_need_pos[i,1],A_need_pos[i,2]]
      out$row[j] <- ifelse(length(rownames(A))>0,rownames(A)[A_need_pos[i,1]],A_need_pos[i,1])
      out$col[j] <- ifelse(length(colnames(A))>0,colnames(A)[A_need_pos[i,2]],A_need_pos[i,2])
      j <- j+1
    }
  }
  if(nrow(S_need_pos)>0){
    for(i in 1:nrow(S_need_pos)){
      Spos[S_need_pos[i,1],S_need_pos[i,2]] <- out$name[j] <- paste(
        model@name,".S[",S_need_pos[i,1],",",S_need_pos[i,2],"]",sep="")
      if(!all.na(model_S$labels)){
        out$label[j] <- model_S$labels[S_need_pos[i,1],S_need_pos[i,2]]
      }
      out$matrix[j] <- "S"
      out$Raw.Value[j] <- S[S_need_pos[i,1],S_need_pos[i,2]]
      out$row[j] <- ifelse(length(rownames(S))>0,rownames(S)[S_need_pos[i,1]],S_need_pos[i,1])
      out$col[j] <- ifelse(length(colnames(S))>0,colnames(S)[S_need_pos[i,2]],S_need_pos[i,2])
      j <- j+1
    }
  }
  #Get standardized values:
  freeparams <- omxGetParameters(model)
  paramnames <- names(freeparams)
  zout <- .standardizeParams(x=freeparams,model=model,Apos=Apos,Spos=Spos)
  #Compute SEs, or assign them 'not requested' values, as the case may be:
  if(SE){ 
    #From Mike Hunter's delta method example:
    covParam <- ParamsCov[paramnames,paramnames]#<--submodel will usually not contain all free param.s
    jacStand <- numDeriv::jacobian(func=.standardizeParams, x=freeparams, model=model, Apos=Apos, Spos=Spos)
    covSparam <- jacStand %*% covParam %*% t(jacStand)
    dimnames(covSparam) <- list(names(zout),names(zout))
    SEs <- sqrt(diag(covSparam))
    #SEs[diag(covSparam)<.Machine$double.eps] <- 0
  }
  else{SEs <- rep("not requested",length(zout)); names(SEs) <- names(zout)}
  #Add standardized values and SEs to output:
  out$Std.Value <- zout
  out$Std.SE <- SEs
  #Pull in raw SEs if requested:
  if(SE){
    for(i in 1:numelem){
      if( (out$name[i] %in% paramnames) | 
            (out$label[i] %in% paramnames) ){
        out$Raw.SE[i] <- sqrt(covParam[ifelse(is.na(out$label[i]),out$name[i],out$label[i]),
                                       ifelse(is.na(out$label[i]),out$name[i],out$label[i])])
  }}}
  else{out$Raw.SE <- "not requested"}
  return(out)
}
mxStandardizeRAMpaths <- function(model, SE=FALSE){
	if (model@.wasRun && model@.modifiedSinceRun){
		msg <- paste("MxModel", omxQuotes(model@name), "was modified",
			     "since it was run.")
		warning(msg)
	}
  #If SE=T,need to check for independent submodels because they will have their own Hessians;
  #recur main function as appropriate:
  inde.subs.flag <- FALSE
  if(SE & length(model$submodels)>0){
    RAM.subs <- (sapply(model@submodels,function(x){class(x$expectation)})=="MxExpectationRAM" | 
                   sapply(model@submodels,function(x){length(x@submodels)>0}))
    inde.subs <- sapply(model@submodels,function(x){x@independent})==TRUE &
      (sapply(model@submodels,function(x){class(x$expectation)})=="MxExpectationRAM" | 
         sapply(model@submodels,function(x){length(x@submodels)>0}))
    if(sum(inde.subs)>0){
      #if ALL submodels are either independent RAM models or non-RAM models:
      if(all(RAM.subs==inde.subs)){ 
        out <- lapply(model@submodels[which(inde.subs)],mxStandardizeRAMpaths,SE=T)
        if(length(out)==0){stop(paste("model '",model@name,"' contains no submodels that use RAM expectation",sep=""))}
        return(out)
      }
      else{out2 <- lapply(model@submodels[which(inde.subs)],mxStandardizeRAMpaths,SE=T)}
    inde.subs.flag <- TRUE
  }}
  covParam <- NULL
  #If user requests SEs, check to be sure they can and should be computed:
  if(SE){
    if(length(model@constraints)>0){
      msg <- paste("standard errors will not be computed because model '",model@name,"' contains at least one mxConstraint",sep="")
      warning(msg)
      SE <- FALSE
    }
    if(SE & length(model@output$hessian)==0){
      if(!model@.wasRun){
        msg <- paste("standard errors will not be computed because model '",model@name,"' has not yet been run",sep="")
        warning(msg)
        SE <- FALSE
      }
      else{
        warning("argument 'SE=TRUE' requires model to have a nonempty 'hessian' output slot; continuing with 'SE' coerced to 'FALSE'")
        SE <- FALSE
    }}
    libraries <- rownames(installed.packages())
    pkgcheck <- ("numDeriv" %in% libraries)
    if(SE & !pkgcheck){
      warning("argument 'SE=TRUE' requires package 'numDeriv' to be installed; continuing with 'SE' coerced to 'FALSE'")
      SE <- FALSE
    }
    if(SE){
      if(model@output$infoDefinite){
        #solve() will fail if Hessian is computationally singular;
        #chol2inv() will still fail if Hessian is exactly singular.
        covParam <- 2*chol2inv(chol(model@output$hessian))
        dimnames(covParam) <- dimnames(model@output$hessian)
      }
      #An indefinite Hessian usually means some SEs will be NaN:
      else{covParam <- 2*solve(model@output$hessian)}
  }}
  #Check if single-group model uses RAM expectation, and proceed if so:
  if(length(model@submodels)==0){
    if(class(model$expectation)!="MxExpectationRAM"){stop(paste("model '",model@name,"' does not use RAM expectation",sep=""))}
    return(.mxStandardizeRAMhelper(model=model,SE=SE,ParamsCov=covParam))
  }
  #Handle multi-group model:
  if(length(model@submodels)>0){
    if(!inde.subs.flag){
      out <- lapply(
        model@submodels[which(
          (sapply(model@submodels,function(x){class(x$expectation)})=="MxExpectationRAM" | 
          sapply(model@submodels,function(x){length(x@submodels)>0}))
          )],
        .mxStandardizeRAMhelper,SE=SE,ParamsCov=covParam)
      if(length(out)==0){stop(paste("model '",model@name,"' does not use RAM expectation",sep=""))}
      return(out)
    }
    else{
      out1 <- lapply(
        model@submodels[which(
          (sapply(model@submodels,function(x){class(x$expectation)})=="MxExpectationRAM" | 
             sapply(model@submodels,function(x){length(x@submodels)>0})) & 
            !sapply(model@submodels,function(x){x@independent})
        )],
        .mxStandardizeRAMhelper,SE=SE,ParamsCov=covParam)
      out <- as.list(c(out1,out2))
      if(length(out)==0){stop(paste("model '",model@name,"' contains no submodels that use RAM expectation",sep=""))}
      out <- out[names(model@submodels[which(
        (sapply(model@submodels,function(x){class(x$expectation)})=="MxExpectationRAM" | 
           sapply(model@submodels,function(x){length(x@submodels)>0}))
      )])]
      return(out)
}}}

coef.MxModel <- function(object, ...) omxGetParameters(object)
