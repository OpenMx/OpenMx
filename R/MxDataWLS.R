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


#------------------------------------------------------------------------------

setClass(Class = "MxDataLegacyWLS",
	 contains = "MxDataStatic",
	 representation = representation(
		acov   = "matrix",
		fullWeight = "matrix",
		thresholds = "matrix"))

setMethod("initialize", "MxDataLegacyWLS",
	  function(.Object, observed, means, type, numObs, acov, fullWeight, thresholds) {
		.Object@observed <- observed
		.Object@means <- means
		.Object@type <- type
		.Object@numObs <- numObs
		.Object@acov <- acov
		.Object@fullWeight <- fullWeight
		.Object@thresholds <- thresholds
		.Object@name <- "data"
		.Object@.needSort <- FALSE
		.Object@.isSorted <- FALSE
		.Object@primaryKey <- as.character(NA)
		.Object@weight <- as.character(NA)
		.Object@frequency <- as.character(NA)
		return(.Object)
	}
)

legacyMxData <- function(observed, type, means = NA, numObs = NA, acov=NA, fullWeight=NA, thresholds=NA) {
	if (length(means) == 1 && is.na(means)) means <- as.numeric(NA)
	if (length(acov) == 1 && is.na(acov)) acov <- matrix(as.numeric(NA))
	if (length(fullWeight) == 1 && is.na(fullWeight)) fullWeight <- matrix(as.numeric(NA))
	if (length(thresholds) == 1 && is.na(thresholds)) thresholds <- matrix(as.numeric(NA))
	if (missing(observed) || !is(observed, "MxDataFrameOrMatrix")) {
		stop("Observed argument is neither a data frame nor a matrix")
	}
	dups <- duplicated(colnames(observed))
	if (any(dups)) {
		stop(paste("Column names must be unique. Duplicated:",
			   omxQuotes(colnames(observed)[dups])))
	}
	if ((!is.vector(means) && !(prod(dim(means)) == length(means))) || !is.numeric(means)) {
		stop("Means argument must be of numeric vector type")
	}
	if (type != "raw" && is.na(numObs)) {
		stop("Number of observations must be specified for non-raw data, i.e., add numObs=XXX to mxData()")
	}
	if (type == "acov") {
		verifyCovarianceMatrix(observed, strictPD=FALSE)
		verifyCovarianceMatrix(acov, nameMatrix="asymptotic", strictPD=FALSE)
		if(!single.na(fullWeight)){
			verifyCovarianceMatrix(fullWeight, nameMatrix="asymptotic", strictPD=FALSE)
		}
		if ( !single.na(thresholds) ) {
			verifyThresholdNames(thresholds, observed)
		}
	}

	lapply(dimnames(observed)[[2]], imxVerifyName, -1)
	if(is.matrix(means)){meanNames <- colnames(means)} else {meanNames <- names(means)}
	means <- as.matrix(means)
	dim(means) <- c(1, length(means))
	colnames(means) <- meanNames

	return(new("MxDataLegacyWLS", observed, means, type, as.numeric(numObs), acov, fullWeight,
		   thresholds))
}

tryCatch.W <- function(expr) {
	# see demo(error.catching)
	W <- NULL
	w.handler <- function(w) {
		W <<- c(W,w)
		invokeRestart("muffleWarning")
	}
	list(value = withCallingHandlers(tryCatch(expr), warning = w.handler),
	     warning = W)
}

#------------------------------------------------------------------------------
# Mike Hunter's wls compute function for continuous only variables

# x is the raw data
wlsContinuousOnlyHelper <- function(x, type="WLS"){
	mnames <- colnames(x)
	numRows <- nrow(x)
	numCols <- ncol(x)
	numColsStar <- numCols*(numCols+1)/2
	if(numRows-1 < numColsStar){
		stop(paste0('Too few rows (', numRows, ') for number of variables (', numCols, ').\nFor WLS, you need at least n*(n+1)/2 + 1 = ', numColsStar+1, ' rows.\nBetter start rubbing two pennies together.'))
	}
	
	if(type=="ULS") {
		useWeight <- diag(1, numColsStar)
	}
	
	x <- x - rep(colMeans(x), each=nrow(x))
	V <- cov(x)*(numRows-1)/numRows
	U <- matrix(0, nrow=numColsStar, ncol=numColsStar)
	
	# Now construct the U MATRIX, from the W ARRAY
	#  Step 1: Generate index matrix M
	row <- 1
	M <- matrix(0, nrow=numColsStar, ncol=2)
	for(j in 1:numCols){
		M[row:(row+numCols-j),] <- cbind(j:numCols, j)
		row <- row+numCols-j+1
	}
	
	if(type=="DLS" || type=="DWLS") {
		for(i in 1:numColsStar){
			U[i,i] <- 1/(sum((x[,M[i, 1]]**2) * (x[,M[i, 2]]**2)) / numRows - V[M[i, 1], M[i, 2]]**2)
		}
		useWeight <- U
	}
	
	#  Step 2: Create the U MATRIX
	for(j in 1:numColsStar){
		for(i in j:numColsStar){
			ind <- c(M[i,], M[j,])
			U[i,j] <- sum(x[,ind[1]] * x[,ind[2]] * x[,ind[3]] * x[,ind[4]]) / numRows - V[ind[1],ind[2]]*V[ind[3],ind[4]]
		}
	}
	U <- vech2full(vech(U))
	fullWeight <- chol2inv(chol(U))
	if(type=="WLS"){
		useWeight <- fullWeight
	}

	nv <- ncol(x)
	covNames <- outer(mnames[1:nv], mnames[1:nv], FUN=paste, sep='_')
	diag(covNames) <- paste0("var_", mnames[1:nv])
	vechs(covNames) <- paste0("poly_", vechs(covNames))
	n1 <- vech(covNames)
	dimnames(useWeight) <- list(n1,n1)
	dimnames(fullWeight) <- list(n1,n1)

	return(list(use=useWeight*numRows, full=fullWeight*numRows))
}




#------------------------------------------------------------------------------
# Ryne Estabrook's wls compute function for only only variables
# Modified by Mike Hunter to allow continuous and maybe joint.


threshLogLik <- function(thresh, rawData, return="individual", useMinusTwo=TRUE){
	# individual: returns -log likelihood for each category 
	#	given particular threshold values
	#	to be used for jacobian
	# model: returns model -log likelihood
	dataTable <- table(rawData)
	thresh <- c(-Inf, thresh, Inf)
	minP   <- pnorm(thresh[1:length(dataTable)])
	maxP   <- pnorm(thresh[2:length(thresh)])
	cellP  <- (- 1 - useMinusTwo) * log(maxP - minP)
	if (return=="individual"){return(cellP)}
	if (return=="model"){return(sum(cellP*dataTable))}
	}

pcLogLik <- function(k, means, vars, thresh, rawData, return="individual", useMinusTwo=TRUE){
	# individual: returns -log likelihood for each category 
	#	given particular threshold values
	#	to be used for jacobian
	# model: returns model -log likelihood
	if (ncol(rawData)!=2)stop("Raw data must contain two variables.")
	if (ncol(thresh)!=2)stop("Threshold matrix must contain two columns.")
	if (length(k)!=1)stop("Please provide a single correlation to be tested.")
	pcThresh  <- rbind(-Inf, thresh, Inf)
	pcThresh[is.na(pcThresh)] <- Inf
	
	# make the frequency counts
	# table() drops missing values
	dataTable <- table(rawData)
	# was 'dataTableTall'
	dtt <- data.frame(
		x=rep(1:dim(dataTable)[1], dim(dataTable)[2]),
		y=rep(1:dim(dataTable)[2], each=dim(dataTable)[1]),
		count=as.vector(dataTable),
		mLL=NA
		)
	dtt$xMin <- pcThresh[dtt$x,	1]
	dtt$xMax <- pcThresh[dtt$x + 1,1]
	dtt$yMin <- pcThresh[dtt$y,    2]
	dtt$yMax <- pcThresh[dtt$y + 1,2]
	
	# make correlation matrix for 
	k <- max(min(k,.999),-.999)
	corMatrix <- matrix(c(1, k, k, 1), 2, 2)
	for (i in 1:dim(dtt)[1]){
		dtt$mLL[i] <- (- 1 - useMinusTwo) * log(mvtnorm::pmvnorm(
			lower=c(dtt$xMin[i], dtt$yMin[i]),
			upper=c(dtt$xMax[i], dtt$yMax[i]),
			mean=c(0, 0),
			corr=corMatrix
			))
		}
	
	if (return=="individual"){return(as.vector(dtt$mLL))}
	if (return=="model"){return(sum(dtt$count*dtt$mLL))}
	if (return=="table"){return(dtt)}
	}

rcLogLik <- function(k, means=NULL, vars=NULL, thresh=NULL, rawData, return="model", useMinusTwo=TRUE){
	if (ncol(rawData)!=2)stop("Raw data must contain exactly two variables.")
	if (length(k)!=1)stop("Please provide a single correlation to be tested.")
	
	if(is.null(means)){ means <- apply(rawData, 2, mean, na.rm=TRUE)}
	if(is.null(vars)){ vars <- apply(rawData, 2, var, na.rm=TRUE)}
	
	sigma <- matrix(c(vars[1], k, k, vars[2]), 2, 2)
	lik <- apply(rawData, 1, mvtnorm::dmvnorm, means, sigma)
	
	if (return=="model"){ return( (- 1 - useMinusTwo)*sum(log(lik), na.rm=TRUE) ) }
	if (return=="individual") { return( (-1-useMinusTwo)*log(lik) ) }
	}

psLogLik <- function(k, means, vars, thresh, rawData, return="model", useMinusTwo=TRUE, print.res=FALSE){
	if (ncol(rawData)!=2)stop("Raw data must contain two variables.")
	if (length(k)!=1)stop("Please provide a single correlation to be tested.")
	isOrd <- unlist(lapply(rawData, is.ordered))
	if (sum(isOrd)!=1)stop("Raw data must contain one ordinal variable and one numeric variable.")
	
	# prep the threshold matrix
	#pcThresh  <- rbind(-Inf, thresh, Inf)
	#pcThresh[is.na(pcThresh)] <- Inf
	
	#wideData <- rawData
	#wideData$min <- pcThresh[as.numeric(wideData[,isOrd])]
	#wideData$max <- pcThresh[as.numeric(wideData[,isOrd])+1]
	
	
	llC <- log(dnorm(rawData[,!isOrd], means[!isOrd], sqrt(vars[!isOrd])))
	#oMean <- (rawData[,!isOrd] - means[!isOrd]) * k / vars[!isOrd]
	#oVar <- vars[isOrd] - k*(1/vars[!isOrd])*k
	z <- ( rawData[,!isOrd] - means[!isOrd] ) / (vars[!isOrd])
	oMean <- k*z
	oVar <- max(c(vars[isOrd] - k*k/vars[!isOrd], 1e-10))
	cumProb <- sapply(thresh, pnorm, mean=oMean, sd=sqrt(oVar))
	cumProb <- matrix(cumProb, nrow=nrow(rawData), ncol=length(thresh))
	cumProb <- cbind(cumProb, 1)
	levProb <- cbind(cumProb[,1,drop=F], cumProb[,-1,drop=F] - cumProb[,-(length(thresh)+1),drop=F])
	sel <- unclass(rawData[,isOrd])
	
	llO <- rep(NA, length(llC))
	for (i in 1:(length(thresh)+1)){
		llO[sel %in% i] <- levProb[sel %in% i, i]
	}
	llO <- log(llO)
	
	if(return=="model"){
		if(print.res) {
			print(paste('k =', k))
			print(paste('-2LL =', (- 1 - useMinusTwo) * sum(llC+llO, na.rm=TRUE)))
		}
		return((- 1 - useMinusTwo) * sum(llC+llO, na.rm=TRUE))
	} else if(return=="individual"){
		return((- 1 - useMinusTwo) * (llC+llO))
	}
	}

normLogLik <- function(pars, rawData, return="model", useMinusTwo=TRUE){
	#ret <- (- 1 - useMinusTwo) * log(dnorm(rawData, pars[1], sqrt(pars[2])))
	ret <- (- 1 - useMinusTwo) * log( 1/(sqrt(2*pi*abs(pars[2]))) * exp(-(rawData-pars[1])^2/(2*abs(pars[2]))))
	if(return=="individual"){
		return(ret)
	}
	if(return=="model"){
		return( sum(ret, na.rm=TRUE) )
	}
}

normLogLikGrad <- function(pars, rawData, return="model", useMinusTwo=TRUE){
	ret <- matrix(NA, nrow=length(rawData), ncol=2)
	ret[,1] <- (rawData - pars[1])/pars[2]
	ret[,2] <- (rawData - pars[1])^2/(2*pars[2]^2) - 1/(2*pars[2])
	ret <- (-1 - useMinusTwo)*ret
	if(return=="individual"){
		return(ret)
	}
	if(return=="model"){
		return( apply(ret, 2, sum) )
	}
}
# Note for the multivariate case
# m = mean parameter vector (not sample mean)
# L = log likelihood
# S = variance parameter matrix
# Sinv = solve(S)
# y = raw data vector for a single row
# r = residual = y - m
#  dL/dm = Sinv %*% (y - m)
#  dL/dS = -0.5*( Sinv - Sinv %*% r %*% t(r) %*% Sinv)


normLogLikHess <- function(pars, rawData, return="model", useMinusTwo=TRUE){
	ret <- matrix(NA, nrow=length(rawData), ncol=2)
	ret[,1] <- (rawData-1)/pars[2]
	ret[,2] <- 1/(2*pars[2]^2) - (rawData-pars[1])^2/(pars[2])^3
	ret <- (-1 - useMinusTwo)*ret
	if(return=="individual"){
		return(ret)
	}
	if(return=="model"){
		return( apply(ret, 2, sum, na.rm=TRUE) )
	}
}
# Note for the multivariate case
# d2L/dm2 = -1*Sinv
# d2L/dS2 = ...
# for 2x2 the Hessian of
# S = a b
#     b c
# is
# 1/(det(S))^2 * 
# c^2       
# -2*b*c    2*(b^2 + a*c)
# b^2       -2*a*b           a^2
#
# D <- matrix(c(1,0,0,0,0,1,1,0,0,0,0,1), 4, 3) #duplication matrix of order 2
# Sinv <- solve(S)
# d2L/dS2 = -0.5 * t(D) %*% ( kronecker(Sinv, Sinv) ) %*% D
# d2L/dSdm = d2L/dmdS = 0
# This is from Abadir and Magnus (2005, p. 390).
#


rc3LogLik <- function(k, means=NULL, vars=NULL, thresh=NULL, rawData, return="model", useMinusTwo=TRUE){
	if (ncol(rawData)!=2)stop("Raw data must contain exactly two variables.")
	if (length(k)!=3)stop("Please provide a variance, a covariance, and another variance to be tested.")
	
	if(is.null(means)){ means <- apply(rawData, 2, mean, na.rm=TRUE)}
	
	sigma <- matrix(c(k[1], k[2], k[2], k[3]), 2, 2)
	lik <- (-1-useMinusTwo)*log(apply(rawData, 1, mvtnorm::dmvnorm, means, sigma))
	
	if (return=="model"){ return(sum(lik)) }
	if (return=="individual") { return(lik) }
}

rc3Hess <- function(k, means=NULL, vars=NULL, thresh=NULL, rawData, return="model", useMinusTwo=TRUE){
	if (ncol(rawData)!=2)stop("Raw data must contain exactly two variables.")
	if (length(k)!=3)stop("Please provide a variance, a covariance, and another variance to be tested.")
	
	if(is.null(means)){ means <- apply(rawData, 2, mean, na.rm=TRUE)}
	
	S <- matrix(c(k[1], k[2], k[2], k[3]), 2, 2)
	D <- matrix(c(1,0,0,0,0,1,1,0,0,0,0,1), 4, 3) #duplication matrix of order 2
	Sinv <- solve2x2(S) #matrix(c(S[2,2], -S[1,2], -S[2,1], S[1,1]), 2, 2)/(S[1,1]*S[2,2] - S[1,2]*S[2,1])
	#Sinv <- solve(S)
	AnalyticCovHessian <- nrow(rawData) * (-1-useMinusTwo) * -0.5 * t(D) %*% ( kronecker(Sinv, Sinv) ) %*% D
	return(AnalyticCovHessian)
}
# Note: 3000x faster (3350x in simulation) than numerical Hessian from numDeriv::hessian of rc3LogLik
# 2000x faster when solve(S) is used instead of solve2x2

indexCov4to2 <- function(i, j, k, l, nvar){
	a <- indexCov2to1(i, j, nvar)
	b <- indexCov2to1(k, l, nvar)
	#indexCov2to1(a, b, nvar*(nvar+1)/2)
	return( c(a, b) ) #return upper triangle element
}

indexCov2to1 <- function(i, j, nvar){
	if(i < j) {stop("Element should be in lower triangle: column j <= row i")}
	a <- i
	b <- j
	return( a + nvar*(b-1) - sum((b-1):0) )
}

solve2x2 <- function(x){
	Xinv <- matrix(c(x[2,2], -x[1,2], -x[2,1], x[1,1]), 2, 2)/(x[1,1]*x[2,2] - x[1,2]*x[2,1])
	return(Xinv)
}

univariateThresholdStatisticsHelper <- function(od, data, nvar, n, ntvar, useMinusTwo){
	### univariate thresholds 
	# prep objects
	nlevel <- unlist(lapply(od, nlevels))
	counts <- lapply(od, table)
	thresh <- matrix(NA, ifelse(nvar > 0, max(nlevel)-1, 0), nvar)
	threshHess <- list(NULL)
	threshWarn <- rep(0, nvar)
	if(nvar > 0) {threshJac <- list(NULL)} else threshJac <- NULL
	
	# get the thresholds, their hessians & their jacobians
	if(nvar > 0){
		for (i in 1:nvar){
			a <- proc.time()
			# threshold & jacobian
			tab <- table(od[,i])
			if(any(tab %in% 0)){
				msg <- paste0("Variable ", omxQuotes(names(od)[i]), " has a zero frequency category ", omxQuotes(names(tab)[tab %in% 0]), ".\nEliminate this level in your mxFactor() or combine categories in some other way.\nDo not pass go. Do not collect $200.")
				stop(msg, call.=FALSE)
			}
			startVals <- qnorm(cumsum(tab)/sum(!is.na(od[,i])))
			if (length(startVals)>2){
				uni <- optim(startVals[1:(length(startVals) - 1)], 
					threshLogLik, return="model", rawData=od[,i], useMinusTwo=useMinusTwo, hessian=TRUE, method="BFGS")
			} else {
				result <- tryCatch.W(optimize(threshLogLik, lower=-6.28, upper=6.28,
							     return="model", rawData=od[,i]))
				threshWarn[i] <- length(result$warning)
				tHold <- result$value
				hHold <- numDeriv::hessian(threshLogLik, x=tHold$minimum, 
					return="model", rawData=od[,i])
				uni <- list(par=tHold$minimum, hessian=hHold)
			}
			# assign thresholds
			thresh[1:(nlevel[i] - 1),i] <- uni$par
			# assign hessians
			threshHess[[i]] <- uni$hessian
			# get jacobian
			jac <- numDeriv::jacobian(func=threshLogLik, x=uni$par, rawData=od[,i])
			# assign jacobian
			threshJac[[i]] <- jac[unclass(od[,i]),]
			proc.time() - a
			}
		threshJac <- matrix(unlist(threshJac), nrow=n)
	}
	names(threshHess) <- names(od)
	colnames(thresh)  <- names(od)
	return(list(thresh, threshHess, threshJac, threshWarn))
}

univariateMeanVarianceStatisticsHelper <- function(ntvar, n, ords, data, useMinusTwo){
	### put the means in!
	
	### Use normLogLik function to get ML estimates of univariate
	# means and variances.
	# And use optim to get Jacobians and Hessians of these ML estimates.
	# Populate correct entries of relevant matrices to be used later.
	# In particular, the pcVars and pcMeans in the bivariate for loops.
	startEst <- numeric(2)
	meanEst <- numeric(ntvar)
	varEst <- numeric(ntvar)
	meanHess <- numeric(ntvar)
	meanJac <- matrix(NA, nrow=n, ncol=ntvar)
	varHess <- numeric(ntvar)
	varJac <- matrix(NA, nrow=n, ncol=ntvar)
	for(i in 1:ntvar){
		if( !ords[i] ){
			startEst[1] <- mean(data[,i], na.rm=TRUE)
			startEst[2] <- var(data[,i], na.rm=TRUE)
			univEst <- optim(par=startEst, fn=normLogLik, rawData=data[,i],
				return="model", useMinusTwo=useMinusTwo,
				method="BFGS", gr=normLogLikGrad, hessian=FALSE)
			# Re-set variance to be positive, if needed.
			if(univEst$par[2] < 0) univEst$par[2] <- -univEst$par[2]
			# N.B. The normal LL function takes the abs() of the variance, so sign flips are possible.
			univHess <- normLogLikHess(pars=univEst$par, rawData=data[,i], return="model", useMinusTwo=useMinusTwo)
			meanEst[i] <- univEst$par[1]
			varEst[i] <- univEst$par[2]
			meanHess[i] <- univHess[1]
			varHess[i] <- univHess[2] #Note: off diagonal elements are analytically exactly zero and are thus discarded.
			univJac <- normLogLikGrad( pars=univEst$par, rawData=data[,i],
					return="individual", useMinusTwo=useMinusTwo) # First column is mean, second is variance
			meanJac[,i] <- univJac[,1]
			varJac[,i] <- univJac[,2]
		} else{
			meanEst[i] <- 0
			varEst[i] <- 1 # Explicitly set ordinal mean and variance to 0 and 1
			meanHess[i] <- 0 #???
			varHess[i] <- 0
			meanJac[,i] <- 0 #???
			varJac[,i] <- rep(0, n)
		}
	}
	return(list(meanEst, varEst, meanHess, varHess, meanJac, varJac))
}

# mxDataWLS itself is deprecated
# useMinusTwo parameter is deprecated
mxDataWLS <- function(data, type="WLS", useMinusTwo=TRUE, returnInverted=TRUE, fullWeight=TRUE,
		      suppressWarnings = TRUE, allContinuousMethod=c("cumulants", "marginals"),
		      silent=!interactive())
{
	allContinuousMethod <- match.arg(allContinuousMethod)
	# version 0.2
	#
	#available types
	wlsTypes <- c("ULS", "DLS", "DWLS", "WLS", "XLS")
	
	# error checking
	if (!is.data.frame(data)){
		stop("'data' must be a data frame.")
	}
	for (cn in colnames(data)) imxVerifyName(cn, 2)
	# check type
	if (!(type %in% wlsTypes)){
		stop(
			paste("Type must be one of '",
				paste(wlsTypes[-length(wlsTypes)], collapse="', '"),
				"', or '", wlsTypes[length(wlsTypes)], "'.", sep="")
			)
	}
	allContinuousMethod <- tolower(allContinuousMethod)
	if (!(allContinuousMethod %in% c('cumulant', 'cumulants', 'marginal', 'marginals'))){
		stop(
			paste("'allContinuousMethod' must be one of ",
				"'cumulants' or 'marginals'.\nBoth plural and singular forms are allowed.", sep="")
			)
	}
	# standardize spelling
	if (allContinuousMethod == 'cumulant') allContinuousMethod <- 'cumulants'
	if (allContinuousMethod == 'marginal') allContinuousMethod <- 'marginals'
		
	# select ordinal variables
	ords <- unlist(lapply(data, is.ordered))
	badFactors <- !ords & unlist(lapply(data, is.factor))
	if (any(badFactors)) {
		stop(paste("Factors", omxQuotes(colnames(data)[badFactors]),
			   "must be ordered and are not"))
	}

	nvar <- sum(ords)
	ntvar <- ncol(data)
	n <- dim(data)[1]
	
	
	msg <- paste("Calculating asymptotic summary statistics for",
		ntvar - nvar, "continuous and", nvar, "ordinal variables ...")
	msgLen <- nchar(msg)
	#message(msg)
	if (!silent) imxReportProgress(msg, 0)
	
	# if no ordinal variables, use continuous-only helper
	if(nvar == 0 && allContinuousMethod %in% c("cumulant", "cumulants")){
		if (!silent) imxReportProgress("", msgLen)
		if (any(is.na(data))) {
			stop(paste("All continuous data with missingness cannot be",
				   "handled in the WLS framework.",
				   "Use na.omit(yourDataFrame) to remove rows with missing values",
				   "or use maximum likelihood instead"))
		}
		wls <- wlsContinuousOnlyHelper(data, type)
		retVal <- legacyMxData(cov(data), type="acov", acov=wls$use, fullWeight=wls$full, numObs=n)
		return(wls.permute(retVal))
	}
	
	# separate ordinal and continuous variables (temporary)
	od <- data[,ords,drop=FALSE]
	cd <- data[,!ords,drop=FALSE]
	
	# thresholds, their hessians & their jacobians
	utsList <- univariateThresholdStatisticsHelper(od, data, nvar, n, ntvar, useMinusTwo)
	thresh <- utsList[[1]]
	threshHess <- utsList[[2]]
	threshJac <- utsList[[3]]
	threshWarn <- utsList[[4]]
	
	# means and variances with their hessians & their jacobians
	umvsList <- univariateMeanVarianceStatisticsHelper(ntvar, n, ords, data, useMinusTwo)
	meanEst <- umvsList[[1]]
	varEst <- umvsList[[2]]
	meanHess <- umvsList[[3]]
	varHess <- umvsList[[4]]
	meanJac <- umvsList[[5]]
	varJac <- umvsList[[6]]
	
	### bivariate polychorics
	pcMatrix <- matrix(NA, ntvar, ntvar)
	diag(pcMatrix) <- 1
	pcJac    <- matrix(NA, nrow=n, ncol=ntvar*(ntvar+1)/2)
	hessHold <- numeric(ntvar*(ntvar+1)/2)
	hessWarn <- rep(0, length(hessHold))
	parName  <- NULL
	r3hess <- array(NA, dim=c(3, 3, ntvar*(ntvar+1)/2))
	covHess <- matrix(0, nrow=ntvar*(ntvar+1)/2, ncol=ntvar*(ntvar+1)/2)
	# det of correlation matrix = 1 - r^2
	# sqrt of det is sqrt(1 - r^2)
	# inverse of correlation matrix is (1-r^2) * [ 1 & -r // -r & 1]
	#   - diagonal elements are 1/(1-r^2)
	#   - off-diagonal elements are -r/(1-r^2)
	# you don't need any of this information now, but may later
	for (j in 1:ntvar){
		for (i in j:ntvar){
			pcData <- data[,c(i,j)]
			ordPair <- (ords[i] + ords[j])
			if( ordPair == 0 ) { # Continuous variables
				logLikFUN <- rcLogLik
				pcThresh <- NULL
			} else if( ordPair == 1 ) { # Joint variables
				logLikFUN <- psLogLik
				# Find correct ordinal column
				tcols <- as.numeric(na.omit(match(names(ords)[c(i,j)], colnames(thresh))))
				pcThresh <- matrix(thresh[,tcols], ncol=1)
				pcThresh <- pcThresh[!is.na(pcThresh),]
			} else if( ordPair == 2 ) { # Ordinal variables
				logLikFUN <- pcLogLik
				# Find correct ordinal column
				tcols <- as.numeric(na.omit(match(names(ords)[c(i,j)], colnames(thresh))))
				pcThresh <- matrix(thresh[,tcols], ncol=2)
			} else stop(paste("Cannot determine variable type for columns", i, "and" , j))
			pcMeans <- c(meanEst[i], meanEst[j])
			pcVars <- c(varEst[i], varEst[j])
			pcBounds <- c(-1, 1) * (sqrt(prod(pcVars)) - 1e-6)
			if ( (i==j) ){
				pcMatrix[i,j] <- pcVars[1]
				pcJac[,indexCov2to1(i,j,ntvar)] <- varJac[,i]
				hessHold[indexCov2to1(i, j,ntvar)] <- varHess[i]
				r3hess[,,indexCov2to1(i, j, ntvar)] <- diag(rep(1, 3))
				covHess[matrix(indexCov4to2(i, j, i, j, ntvar), ncol=2)] <- varHess[i]
				parName <- c(parName, paste("var", names(pcData)[1], sep="_"))
			} else {
				# parameter name
				parName <- c(parName, paste("poly", names(pcData)[1], names(pcData)[2], sep="_"))
				# get polychoric
				optResult <- tryCatch.W(optimize(logLikFUN, lower=pcBounds[1], upper=pcBounds[2],
							means=pcMeans, vars=pcVars, thresh=pcThresh, return="model", rawData=pcData, useMinusTwo=useMinusTwo))
				pc <- optResult$value
				hessWarn[indexCov2to1(i, j, ntvar)] <- length(optResult$warning)
				# assign polychoric
				pcMatrix[j, i] <- pc$minimum
				pcMatrix[i, j] <- pc$minimum
				# get and assign hessian
				# Stop Hessian from walking outside of bounds
				small <- 0.1
				step <- pc$minimum + c(-1, 1)*.1
				if(pcBounds[1] > step[1] || pcBounds[2] < step[2]){
					# if we're within 0.1 of the bound then only walk halfway to it.
					small <- min(abs(pcBounds - pc$minimum))/2
				}
				# Compute actual Hessian
				hessHold[indexCov2to1(i, j, ntvar)] <- numDeriv::hessian(logLikFUN, x=pc$minimum, 
						means=pcMeans, vars=pcVars, thresh=pcThresh, return="model", rawData=pcData, useMinusTwo=useMinusTwo, method.args=list(d=small))
#				if(ordPair==0){ #Continuous variables
#					r3hess[,,indexCov2to1(i, j, ntvar)] <- numDeriv::hessian(rc3LogLik, x=c(pcVars[1], pc$minimum, pcVars[2]),
#						means=meanEst[c(i, j)], thresh=pcThresh, return="model", rawData=pcData, useMinusTwo=useMinusTwo)
#					r3hess[,,indexCov2to1(i, j, ntvar)] <- cov2cor(r3hess[,,indexCov2to1(i, j, ntvar)])
#					covHess[matrix(indexCov4to2(i, i, i, j, ntvar), ncol=2)] <- r3hess[1,2,indexCov2to1(i, j, ntvar)]*sqrt(varHess[i]*hessHold[indexCov2to1(i, j, ntvar)])
#					#covHess[matrix(indexCov4to2(i, i, j, j, ntvar), ncol=2)] <- r3hess[1,3,indexCov2to1(i, j, ntvar)]*sqrt(varHess[i]*varHess[j])
#					covHess[matrix(indexCov4to2(i, j, j, j, ntvar), ncol=2)] <- r3hess[2,3,indexCov2to1(i, j, ntvar)]*sqrt(hessHold[indexCov2to1(i, j, ntvar)]*varHess[j])
#				}
				#pcHess[i, j]   <- pcHess[j, i]
				# get jacobian
				if( ordPair == 0 ) { # Continuous variables
					assignJac <- matrix(numDeriv::jacobian(func=logLikFUN, x=pc$minimum, means=pcMeans, vars=pcVars, thresh=pcThresh, 
						rawData=pcData, return="individual", useMinusTwo=useMinusTwo), 
						nrow=n, 
						ncol=1)
				}
				else if( ordPair == 1 ) { #Joint variables
					assignJac <- matrix(numDeriv::jacobian(func=logLikFUN, method.args=list(eps=1e-4, d=1e-4, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE), x=pc$minimum, means=pcMeans, vars=pcVars, thresh=pcThresh, 
						rawData=pcData, return="individual", useMinusTwo=useMinusTwo), 
						nrow=n, 
						ncol=1)
					#stop("Jacobian for joint ordinal and continuous variables is not yet implemented.")
				} else { # Ordinal variables
					localJac <- matrix(numDeriv::jacobian(func=logLikFUN, x=pc$minimum, vars=pcVars, thresh=pcThresh, 
						rawData=pcData, return="individual", useMinusTwo=useMinusTwo), 
						nrow=nlevels(pcData[,1]), 
						ncol=nlevels(pcData[,2]))
					# For all ordinal data this is the n1 x n2 table of the 1x1 gradients in each combination
					#  of levels of the ordered variables.
					# Then the following loop turns the matrix into a vector corresponding to the rows of data.
					# The continuous data assignJac already returns the vector of 1x1 gradients
					#  (correponding to data rows).
					# assign jacobian
					select <- matrix(c(unclass(pcData[,1]), unclass(pcData[, 2])), ncol=2)
					assignJac <- localJac[select]
					# N.B. when either variable is missing, this returns NA for that row.
					}
					pcJac[,indexCov2to1(i,j,ntvar)] <- assignJac
				}
			}
		}
	# nparam <- nvar*(nvar+1)/2 + length(threshHess) #N.B. not used anywhere
	
	
	covHess[lower.tri(covHess)] <- t(covHess)[lower.tri(covHess)]
	diag(covHess) <- hessHold
	colnames(covHess) <- parName
	rownames(covHess) <- parName
	
	### put names on everything
	
	# put names on the polychorics
	dimnames(pcMatrix) <- list(names(data), names(data))
	
	# put names on the jacobian
	colnames(pcJac) <- parName
	
	# put names on the hessian
	names(hessHold) <- parName
	
	# Now doing something about the means!
	# To replicate old behavior set,
	# The following two lines should be deleted.
	#meanJac <- NULL
	#meanHess <- NULL
	# even though these might not be NULL and have been processed earlier.
	fullJac  <- cbind(pcJac, meanJac, threshJac)
	if( nvar > 0 ){
		fullHess <- as.matrix(Matrix::bdiag(diag(c(hessHold, meanHess)), Matrix::bdiag(threshHess)))
	} else {
		fullHess <- diag(c(hessHold, meanHess))
	}
	
	#TODO Figure out why certain elements of fullJac end up missing when the data are missing.
	quad <- (n-1)*var(fullJac, use="pairwise.complete.obs") #bc colMeans all zero == t(fullJac) %*% fullJac
	quad[is.na(quad)] <- 0
	sel  <- diag(quad)!=0
	iqj  <- matrix(0, dim(quad)[1], dim(quad)[2])
	attIqj <- try(solve(quad[sel, sel]))
	if(class(attIqj) %in% "try-error"){
		iqj[sel,sel] <- MASS::ginv(quad[sel, sel])
		warning('First derivative matrix was not intertible. Used pseudo-inverse instead.')
	} else {
		iqj[sel,sel] <- attIqj
	}
	
	# make the weight matrix!!!
	wls <- fullHess %*% iqj %*% fullHess
	fullNames <- c(parName, colnames(data))
	if(nvar > 0){
		fullNames <- c(fullNames,
				unlist(mapply(function(vn,hess) paste0(vn,'t',1:ncol(hess)),
					names(threshHess),
					threshHess)))
	}
	dimnames(wls) <- list(fullNames, fullNames)
	fullWarn <- c(hessWarn, threshWarn)
	names(fullWarn) <- c(parName, colnames(data)[ords])
	if (!suppressWarnings && any(fullWarn > 0)) {
		warning(paste("Encountered warnings during optimization of",
			      omxQuotes(fullWarn[fullWarn > 0])))
	}

	dls <- diag(diag(wls))
	dimnames(dls) <- dimnames(wls)
	uls <- (dls>0)*1
	dimnames(uls) <- dimnames(wls)

	# try the weird non-hao version
	xls <- quad
	dimnames(xls) <- dimnames(wls)
	
	dummy <- diag(1, nrow=nrow(pcMatrix))
	dimnames(dummy) <- dimnames(pcMatrix)
	if(nvar > 0){
		retVal <- legacyMxData(dummy, type="acov", numObs=n, 
			acov=diag(1), fullWeight=NA, thresholds=thresh)
	} else {
		retVal <- legacyMxData(dummy, type="acov", numObs=n, 
			acov=diag(1), fullWeight=NA, thresholds=NA)
	}
	retVal@observed <- pcMatrix
	if(fullWeight==TRUE){
		retVal@fullWeight <- wls
	}
	retVal@means <- matrix(meanEst, nrow=1)
	dimnames(retVal@means) <- list(NULL, names(data))
	if (type=="ULS"){
		retVal@acov <- uls
		}
	if (type=="DLS" || type=="DWLS"){
		retVal@acov <- dls
		}	
	if (type=="WLS"){
		retVal@acov <- wls
		}
	if (type=="XLS"){
		retVal@acov <- xls
		}
	if (!silent) imxReportProgress("", msgLen)
	return(wls.permute(retVal))
}

wls.permute <- function(mxd) {
	acov <- mxd@acov
	perm <- match(names(.mxDataAsVector(mxd)), colnames(acov))
	mxd@acov <- acov[perm,perm]
	fw <- mxd@fullWeight
	if (!is.null(fw)) {
		mxd@fullWeight <- fw[perm,perm]
	}
	mxd
}

.mxDataAsVector <- function(mxd) {
	mnames <- colnames(mxd@observed)
	ordInd <- c()
	if (!is.null(mxd@thresholds)) {
		ordInd <- match(colnames(mxd@thresholds), mnames)
		dth <- !is.na(mxd@thresholds)
	}
	v <- c()
	vn <- c()
	for (vx in 1:length(mnames)) {
		tcol <- which(vx == ordInd)
		if (length(tcol) == 0) {
			if (!single.na(mxd@means)) {
				v <- c(v, mxd@means[vx])
				vn <- c(vn, mnames[vx])
			}
		} else {
			tcount <- sum(dth[,tcol])
			v <- c(v, mxd@thresholds[1:tcount,tcol])
			vn <- c(vn, paste0(mnames[vx], 't', 1:tcount))
		}
	}
	for (vx in 1:length(mnames)) {
		if (any(vx == ordInd)) next
		v <- c(v, mxd@observed[vx,vx])
		vn <- c(vn, paste0('var_', mnames[vx]))
	}
	v <- c(v, vechs(mxd@observed))
	nv <- length(mnames)
	vn <- c(vn, paste0('poly_', vechs(outer(mnames[1:nv], mnames[1:nv], FUN=paste, sep='_'))))
	names(v) <- vn
	v
}


##' Estimate summary statistics used by the WLS fit function
##'
##' The summary statistics are returned in the observedStats slot of
##' the MxData object.
##'
##' @param mxd an MxData object containing raw data
##' @param type the type of WLS weight matrix
##' @param allContinuousMethod which method to use when all indicators are continuous
##' @param ...  Not used.  Forces remaining arguments to be specified by name.
##' @param exogenous names variables to be modelled as exogenous
##' @param fullWeight whether to produce a fullWeight matrix
##' @param returnModel whether to return the whole mxModel (TRUE) or just the mxData (FALSE)
##' @param silent logical. Whether to print status to terminal.
##' @seealso
##' \link{mxFitFunctionWLS}
##' @examples
##' omxAugmentDataWithWLSSummary(mxData(Bollen[,1:8], 'raw'))

omxAugmentDataWithWLSSummary <- function(mxd, type=c('WLS','DWLS','ULS'),
			      allContinuousMethod=c("cumulants", "marginals"),
			      ..., exogenous=c(), fullWeight=TRUE, returnModel=FALSE,
			      silent=TRUE)
{
	type <- match.arg(type)
	allContinuousMethod <- match.arg(allContinuousMethod)
  prohibitDotdotdot(list(...))
	if (mxd@type != 'raw') stop("Data must contain a raw data frame")
	data <- mxd@observed
	notFound <- is.na(match(exogenous, colnames(data)))
	if (any(notFound)) {
		stop(paste("Cannot find exogenous", omxQuotes(exogenous), "in data"))
	}
	nc <- ncol(data)
	names <- setdiff(colnames(data), exogenous)
	weight <- mxd@weight
	if (!is.na(weight)) {
		nc <- nc - 1
		names <- names[-match(weight, colnames(data))]
	}
	frequency <- mxd@frequency
	if (!is.na(frequency)) {
		nc <- nc - 1
		names <- names[-match(frequency, colnames(data))]
	}
	numManifests <- nc - length(exogenous)
	fake <- mxModel("fake",
		mxd,
		mxMatrix("Full", nc, nc, dimnames=list(c(names,exogenous),
			c(names,exogenous)), name="S"),
		mxMatrix("Full", nc, nc, name="A"),
		mxMatrix("Full", numManifests, nc,
			dimnames=list(names,c(names,exogenous)), name="F"),
		mxExpectationRAM(),
		mxFitFunctionWLS(type, allContinuousMethod, fullWeight),
		mxComputeOnce('fitfunction', 'fit'))

	if (length(exogenous)) fake$A$values[1:numManifests,(numManifests+1):nc] <- 1
	fake$S$values[1:numManifests,1:numManifests] <- diag(numManifests)
	fake$F$values[1:numManifests,1:numManifests] <- diag(numManifests)

	ords <- unlist(lapply(data, is.ordered)) & colnames(data) %in% names

	if (any(ords) || allContinuousMethod != 'cumulants') {
		fake <- mxModel(fake,
			mxMatrix(values=0, nrow=1, ncol=nc,
				dimnames=list(c(), c(names,exogenous)), name="M"))
		if (length(exogenous)) {
			fake$M$labels[1,(numManifests+1):ncol(fake$M)] <-
				paste0('data.',exogenous)
		}
		fake$expectation$M <- "M"
	}

	if (any(ords)) {
		nthr <- sapply(data[,ords], nlevels) - 1L
		tmpThr <- matrix(NA, ncol=sum(ords), nrow=max(nthr))
		colnames(tmpThr) <- colnames(data)[ords]
		for (cx in 1:ncol(tmpThr)) {
			tmpThr[1:nthr[cx],cx] <- seq(-1,1,length.out=nthr[cx])
		}
		fake <- mxModel(fake, mxMatrix(values=tmpThr, name="thresh"))
		fake$expectation$thresholds <- "thresh"
	}
	fake <- mxRun(fake, silent=silent)
	if (returnModel) return(fake)
	fake$data
}

dataIsRawish <- function(mxd) mxd@type == 'raw' || mxd@type == 'acov'
