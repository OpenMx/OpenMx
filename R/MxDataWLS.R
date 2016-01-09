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


#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Mike Hunter's wls compute function for continuous only variables

# x is the raw data
wlsContinuousOnlyHelper <- function(x, type="WLS"){
	numRows <- nrow(x)
	numCols <- ncol(x)
	numColsStar <- numCols*(numCols+1)/2
	
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
	
	if(type=="DWLS") {
		for(i in 1:numColsStar){
			U[i,i] <- 1/(sum((x[,M[i, 1]]**2) * (x[,M[i, 2]]**2)) / numRows - V[M[i, 1], M[i, 2]]**2)
		}
		useWeight <- U
	}
	
	#Step 2: Create the U MATRIX
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
	
	if (return=="model"){ return( (- 1 - useMinusTwo)*sum(log(lik)) ) }
	if (return=="individual") { return( (-1-useMinusTwo)*log(lik) ) }
	}

psLogLik <- function(k, means, vars, thresh, rawData, return="model", useMinusTwo=TRUE){
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
	z <- ( rawData[,!isOrd] - means[!isOrd] ) / sqrt(vars[!isOrd])
	oMean <- k*z
	oVar <- vars[isOrd] - k*k
	cumProb <- sapply(thresh, pnorm, mean=oMean, sd=sqrt(oVar))
	cumProb <- cbind(cumProb, 1)
	levProb <- cbind(cumProb[,1], cumProb[,-1] - cumProb[,-(length(thresh)+1)])
	sel <- as.numeric(rawData[,isOrd])
	
	llO <- rep(NA, length(llC))
	for (i in 1:(length(thresh)+1)){
		llO[sel==i] <- levProb[sel==i, i]
	}
	llO <- log(llO)
	
	if(return=="model"){
		return((- 1 - useMinusTwo) * sum(llC+llO))
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
		return( sum(ret) )
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
		return( apply(ret, 2, sum) )
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
	if(i > j) {stop("Element should be in upper triangle: column j >= row i")}
	a <- j
	b <- i
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
	if(nvar > 0){
		varn <- apply(!is.na(od), 2, sum)
	}
	counts <- lapply(od, table)
	thresh <- matrix(NA, ifelse(nvar > 0, max(nlevel)-1, 0), nvar)
	threshHess <- list(NULL)
	if(nvar > 0) {threshJac <- list(NULL)} else threshJac <- NULL
	
	# get the thresholds, their hessians & their jacobians
	if(nvar > 0){
		for (i in 1:nvar){
			a <- proc.time()
			# threshold & jacobian
			startVals <- qnorm(cumsum(table(data[,i]))/sum(!is.na(data[,i])))
			if (length(startVals)>2){
				uni <- optim(startVals[1:(length(startVals) - 1)], 
					threshLogLik, return="model", rawData=data[,i], useMinusTwo=useMinusTwo, hessian=TRUE, method="BFGS")
				} else {
				tHold <- optimize(threshLogLik, lower=-6.28, upper=6.28,
					return="model", rawData=data[,i])
				hHold <- numDeriv::hessian(threshLogLik, x=tHold$minimum, 
					return="model", rawData=data[,i])
				uni <- list(par=tHold$minimum, hessian=hHold)
				}
			# assign thresholds
			thresh[1:(nlevel[i] - 1),i] <- uni$par
			# assign hessians
			threshHess[[i]] <- uni$hessian
			# get jacobian
			jac <- numDeriv::jacobian(func=threshLogLik, x=uni$par, rawData=data[,i])
			# assign jacobian
			threshJac[[i]] <- jac[unclass(data[,i]),]
			proc.time() - a
			}
		threshJac <- matrix(unlist(threshJac), nrow=n)
	}
	names(threshHess) <- names(od)
	colnames(thresh)  <- names(od)
	return(list(thresh, threshHess, threshJac))
}

univariateMeanVarianceStatisticsHelper <- function(ntvar, n, ords, data, useMinusTwo){
	### put the means in!
	# means are missing for ordinal data, so 
	meanHess <- NULL
	meanJac <- NULL
	
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

mxDataWLS <- function(data, type="WLS", useMinusTwo=TRUE, returnInverted=TRUE, debug=FALSE, fullWeight=TRUE){
	message("Calculating asymptotic summary statistics ...")
	# version 0.2
	#
	#available types
	wlsTypes <- c("ULS", "DLS", "WLS", "XLS")
	
	# error checking
	if (!is.data.frame(data)){
		stop("'data' must be a data frame.")
		}
	# check type
	if (!(type %in% wlsTypes)){
		stop(
			paste("Type must be one of '",
				paste(wlsTypes[-length(wlsTypes)], collapse="', '"),
				"', or '", wlsTypes[length(wlsTypes)], "'.", sep="")
			)
		}
		
	# select ordinal variables
	ords <- unlist(lapply(data, is.ordered))
	nvar <- sum(ords)
	ntvar <- ncol(data)
	n <- dim(data)[1]

	# if no ordinal variables, use continuous-only helper
	if(nvar ==0){ #N.B. This fails for any missing data
		wls <- wlsContinuousOnlyHelper(data, type)
		return(mxData(cov(data), type="acov", acov=wls$use, fullWeight=wls$full, numObs=n))
	}
	
	# separate ordinal and continuous variables (temporary)
	od <- data[,ords]
	cd <- data[,!ords]
	
	if(ncol(od) > 0 && ncol(cd) > 0 && type=="WLS"){
		msg <- paste("Both ordinal and coninuous variables found with type='WLS'.\n",
			"There is currently a bug that produces biased estimates in this case.\n",
			"Try using type='ULS' or type='DLS'. There is no bug there.\n")
		warning(msg, call.=FALSE)
	}
	
	# thresholds, their hessians & their jacobians
	utsList <- univariateThresholdStatisticsHelper(od, data, nvar, n, ntvar, useMinusTwo)
	thresh <- utsList[[1]]
	threshHess <- utsList[[2]]
	threshJac <- utsList[[3]]
	
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
	parName  <- NULL
	r3hess <- array(NA, dim=c(3, 3, ntvar*(ntvar+1)/2))
	covHess <- matrix(0, nrow=ntvar*(ntvar+1)/2, ncol=ntvar*(ntvar+1)/2)
	# det of correlation matrix = 1 - r^2
	# sqrt of det is sqrt(1 - r^2)
	# inverse of correlation matrix is (1-r^2) * [ 1 & -r // -r & 1]
	#   - diagonal elements are 1/(1-r^2)
	#   - off-diagonal elements are -r/(1-r^2)
	# you don't need any of this information now, but may later
	for (i in 1:ntvar){
		for (j in i:ntvar){
			pcData <- data[,c(i,j)]
			ordPair <- (ords[i] + ords[j])
			if( ordPair == 0 ) { # Continuous variables
				logLikFUN <- rcLogLik
				pcThresh <- NULL
			} else if( ordPair == 1 ) { # Joint variables
				logLikFUN <- psLogLik
				pcThresh <- matrix(thresh[,ifelse(ords[i], i, j)], ncol=1)
				pcThresh <- pcThresh[!is.na(pcThresh),]
			} else if( ordPair == 2 ) { # Ordinal variables
				logLikFUN <- pcLogLik
				pcThresh <- matrix(thresh[,c(i,j)], ncol=2)
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
				pc <- optimize(logLikFUN, lower=pcBounds[1], upper=pcBounds[2],
					means=pcMeans, vars=pcVars, thresh=pcThresh, return="model", rawData=pcData, useMinusTwo=useMinusTwo)
				# assign polychoric
				pcMatrix[j, i] <- pc$minimum
				pcMatrix[i, j] <- pc$minimum
				# get and assign hessian
				hessHold[indexCov2to1(i, j, ntvar)] <- numDeriv::hessian(logLikFUN, x=pc$minimum, 
						means=pcMeans, vars=pcVars, thresh=pcThresh, return="model", rawData=pcData, useMinusTwo=useMinusTwo)
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
	
	# put it all together
	if(nvar == 0){ #jumps to continuousonly ADF code.
		w3 <- diag(hessHold) %*% solve(t(pcJac)%*%pcJac) %*% diag(hessHold)
		w2 <- covHess %*% solve(t(pcJac)%*%pcJac) %*% covHess
		w <- wlsContinuousOnlyHelper(cd)
	}
	# TODO: something still must be done about the means!
	# To replicate old behavior set
	# The following two lines should be deleted.
	meanJac <- NULL
	meanHess <- NULL
	# even though these might not be NULL and have been processed earlier.
	fullJac  <- cbind(pcJac, meanJac, threshJac)
	if( nvar > 0 ){
		fullHess <- as.matrix(Matrix::bdiag(diag(c(hessHold, meanHess)), Matrix::bdiag(threshHess)))
	} else {
		fullHess <- diag(c(hessHold, meanHess))
	}
	
	#TODO Figure out why certain elements of fullJac end up missing when the data are missing.
	quad <- (n-1)*var(fullJac, use="pairwise.complete.obs") #bc colMeans all zero == t(fullJac) %*% fullJac
	sel  <- diag(quad)!=0
	iqj  <- matrix(0, dim(quad)[1], dim(quad)[2])
	iqj[sel,sel] <- solve(quad[sel, sel])
	
	# make the weight matrix!!!
	wls <- fullHess %*% iqj %*% fullHess
	dls <- diag(diag(wls))
	uls <- (dls>0)*1
	
	# try the weird non-hao version
	xls <- quad
	
	if(debug){
		custom.compute <- mxComputeSequence(list(mxComputeNumericDeriv(checkGradient=FALSE), mxComputeReportDeriv()))
		satModel <- mxModel('SatModel4Hess',
			mxMatrix('Symm', nrow(pcMatrix), ncol(pcMatrix), values=pcMatrix, free=TRUE, name='theCov'),
			mxMatrix('Full', 1, length(meanEst), values=meanEst, free=TRUE, name='theMeans'),
			mxMatrix('Full', nrow(thresh), ncol(thresh), values=thresh, free=TRUE, name='theThresh'),
			mxExpectationNormal(covariance='theCov', means='theMeans', thresholds='theThresh', dimnames=names(data)),
			mxFitFunctionML(),
			mxData(data, 'raw'),
			custom.compute)
		
		run <- mxRun(satModel)
		mlHess <- run$output$hessian
		meanID <- grep('.theMeans', rownames(mlHess))
		
		retVal2 <- mxData(pcMatrix, type="acov", numObs=n, 
			acov=diag(1), fullWeight=NA, thresholds=thresh)
		retVal2@acov <- satModel$output$hessian
	}
	
	if(fullWeight==TRUE){
		fw <- wls
	} else {fw <- NA}
	retVal <- mxData(pcMatrix, type="acov", numObs=n, 
		acov=diag(1), fullWeight=NA, thresholds=thresh)
	retVal@fullWeight <- fw
	if (type=="ULS"){
		retVal@acov <- uls
		}
	if (type=="DLS"){
		retVal@acov <- dls
		}	
	if (type=="WLS"){
		retVal@acov <- wls
		}
	if (type=="XLS"){
		retVal@acov <- xls
		}
	if (debug){return(list(fullJac, fullHess))}	
	return(retVal)
	}


