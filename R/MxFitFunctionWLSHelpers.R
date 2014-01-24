#
#   Copyright 2007-2014 The OpenMx Project
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
omxWeightMatrix <- function(type, data, use){
	print("This function is not yet implemented.")
}


#------------------------------------------------------------------------------
# Mike Hunter's wls compute function for continuous only variables

# x is the raw data
wlsContinuousOnlyHelper <- function(x, type="WLS"){
	numRows <- nrow(x)
	numCols <- ncol(x)
	numColsStar <- numCols*(numCols+1)/2
	
	if(type=="ULS") {return(diag(1, numColsStar))}
	
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
		return(U)
	}
	
	#Step 2: Create the U MATRIX
	for(j in 1:numColsStar){
		for(i in j:numColsStar){
			ind <- c(M[i,], M[j,])
			U[i,j] <- sum(x[,ind[1]] * x[,ind[2]] * x[,ind[3]] * x[,ind[4]]) / numRows - V[ind[1],ind[2]]*V[ind[3],ind[4]]
		}
	}
	U <- vech2full(vech(U))
	return(chol2inv(chol(U)))
}




#------------------------------------------------------------------------------
# Ryne Estabrook's wls compute function for only only variables


threshLogLik <- function(thresh, rawData, return="individual", useMinusTwo=FALSE){
	# individual: returns -log likelihood for each category 
	#    given particular threshold values
	#    to be used for jacobian
	# model: returns model -log likelihood
	dataTable <- table(rawData)
	thresh <- c(-Inf, thresh, Inf)
	minP   <- pnorm(thresh[1:length(dataTable)])
	maxP   <- pnorm(thresh[2:length(thresh)])
	cellP  <- (- 1 - useMinusTwo) * log(maxP - minP)
	if (return=="individual"){return(cellP)}
	if (return=="model"){return(sum(cellP*dataTable))}
	}
	
pcLogLik <- function(k, thresh, rawData, return="individual", useMinusTwo=FALSE){
	# individual: returns -log likelihood for each category 
	#    given particular threshold values
	#    to be used for jacobian
	# model: returns model -log likelihood
	require(mvtnorm)
	if (ncol(rawData)!=2)stop("Raw data must contain two variables.")
	if (ncol(thresh)!=2)stop("Threshold matrix must contain two columns.")
	if (length(k)!=1)stop("Please provide a single correlation to be tested.")
	pcThresh  <- rbind(-Inf, thresh, Inf)
	pcThresh[is.na(pcThresh)] <- Inf
	
	# make the frequency counts
	dataTable <- table(rawData)
	# was 'dataTableTall'
	dtt       <- data.frame(
		x=rep(1:dim(dataTable)[1], dim(dataTable)[2]),
		y=rep(1:dim(dataTable)[2], each=dim(dataTable)[1]),
		count=as.vector(dataTable),
		mLL=NA
		)
	dtt$xMin <- pcThresh[dtt$x,    1]
	dtt$xMax <- pcThresh[dtt$x + 1,1]
	dtt$yMin <- pcThresh[dtt$y,    2]
	dtt$yMax <- pcThresh[dtt$y + 1,2]
		
	# make correlation matrix for 
	corMatrix <- matrix(c(1, k, k, 1), 2, 2)
	for (i in 1:dim(dtt)[1]){
		dtt$mLL[i] <- (- 1 - useMinusTwo) * log(pmvnorm(
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

wlsOrdinalOnlyHelper <- function(data, type="WLS", useMinusTwo=FALSE){
	# version 0.2
	#
	
	# load required libraries
	require(OpenMx)
	require(mvtnorm)
	require(numDeriv)
	
	# error checking
	if (!is.data.frame(data)){
		stop("'data' must be a data frame with ordinal variables (factors).")
		}
	# check type
	# must be a better way
	# MDH: Put in a better way
	if (!(type %in% c("ULS", "DLS", "WLS", "XLS"))){
		stop("Type must be either 'ULS', 'DLS', 'WLS', 'XLS'")
		}
		
	# select ordinal variables
	ords <- unlist(lapply(data, is.ordered))
	nvar <- sum(ords)
	
	# discard non-ordinal variables (temporary)
	od      <- data[,ords]
	
	### univariate thresholds 
	# prep objects
	nlevel <- unlist(lapply(od, nlevels)) 
	varn   <- apply(!is.na(od), 2, sum)
	counts <- lapply(od, table)
	thresh <- matrix(NA, max(nlevel)-1, nvar)
	threshHess       <- list(NULL)
	threshJac        <- list(NULL)
	
	# get the thresholds, their hessians & their jacobians
	for (i in 1:nvar){
		a <- proc.time()
		# threshold & jacobian
		startVals <- qnorm(cumsum(table(data[,i]))/sum(!is.na(data[,i])))
		print(startVals)
		print(i)
		if (length(startVals)>2){
			uni <- optim(startVals[1:(length(startVals) - 1)], 
				threshLogLik, return="model", rawData=data[,i], useMinusTwo=useMinusTwo, hessian=TRUE)
			} else {
			tHold <- optimize(threshLogLik, lower=-6.28, upper=6.28,
				return="model", rawData=data[,i])
			hHold <- hessian(threshLogLik, x=tHold$minimum, 
				return="model", rawData=data[,i])
			uni <- list(par=tHold$minimum, hessian=hHold)
			}
		# assign thresholds
		thresh[1:(nlevel[i] - 1),i] <- uni$par
		# assign hessians
		threshHess[[i]] <- uni$hessian
		# get jacobian
		jac <- jacobian(func=threshLogLik, x=uni$par, rawData=data[,i])
		# assign jacobian
		threshJac[[i]] <- jac[unclass(data[,i]),]
		proc.time() - a
		}
	names(threshHess) <- names(od)
	names(threshJac)  <- names(od)
	colnames(thresh)  <- names(od)
	
	### bivariate polychorics
	pcMatrix <- matrix(NA, nvar, nvar)
	diag(pcMatrix) <- 1
	nparam   <- sum(lower.tri(pcMatrix))
	pcJac    <- NULL
	hessHold <- NULL
	parName  <- NULL
	# det of correlation matrix = 1 - r^2
	# sqrt of det is sqrt(1 - r^2)
	# inverse of correlation matrix is (1-r^2) * [ 1 & -r // -r & 1]
	#   - diagonal elements are 1/(1-r^2)
	#   - off-diagonal elements are -1/(1-r^2)
	# you don't need any of this information now, but may later
	
	for (i in 1:(nvar-1)){
		for (j in (i+1):nvar){
			# data
			pcData <- data[,c(i,j)]
			pcThresh <- thresh[,c(i,j)]
			# parameter name
			parName <- c(parName, paste("poly", names(pcData)[1], names(pcData)[2], sep="_"))
			# get polychoric
			pc <- optimize(pcLogLik, lower=-.99, upper=.99, 
				thresh=pcThresh, return="model", rawData=pcData, useMinusTwo=useMinusTwo)
			# assign polychoric
			pcMatrix[j, i] <- pc$minimum
			pcMatrix[i, j] <- pc$minimum
			# get and assign hessian
			hessHold <- c(hessHold, 
				hessian(pcLogLik, x=pc$minimum, 
				thresh=pcThresh, return="model", rawData=pcData, useMinusTwo=useMinusTwo))
			#pcHess[i, j]   <- pcHess[j, i]
			# get jacobian
			localJac <- matrix(jacobian(func=pcLogLik, x=pc$minimum, thresh=pcThresh, 
				rawData=pcData, return="individual"), 
				nrow=nlevels(pcData[,1]), 
				ncol=nlevels(pcData[,2]))
			# assign jacobian
			assignJac <- rep(NA, dim(data)[1])
			cleanData <- unclass(pcData)
			for (k in 1:nlevels(pcData[,2])){
				select <- unclass(pcData[,2])==k
				assignJac[select] <- localJac[,k][unclass(pcData[select,1])]
				}
			pcJac <- cbind(pcJac, assignJac)
			}
		}
	
	# make the hessian
	pcHessI <- diag(nparam)
	invHessULS <- pcHessI
	diag(pcHessI) <- 1/hessHold
	invHessDLS <- pcHessI
	
	### put names on everything
	
	# put names on the polychorics
	dimnames(pcMatrix) <- list(names(od), names(od))
	
	# put names on the jacobian
	colnames(pcJac) <- parName
	
	# put names on the hessian
	dimnames(pcHessI) <- list(parName, parName)
	
	# make the weight matrix!!!
	invHessWLS <- (pcHessI) %*% t(pcJac) %*% pcJac %*% (pcHessI)
	
	if (type=="ULS"){
		return(list(polychoric=pcMatrix, thresh=thresh, invHess=invHessULS))
		}
	if (type=="DLS"){
		return(list(polychoric=pcMatrix, thresh=thresh, invHess=invHessDLS))
		}	
	if (type=="WLS"){
		return(list(polychoric=pcMatrix, thresh=thresh, invHess=invHessWLS))
		}
	}







