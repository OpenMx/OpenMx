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


#--------------------------------------------------------------------
# Author: Michael D. Hunter
# Filename: MxExpectationStateSpace.R
# Date: 2012.11.14
# Purpose: Define classes and methods for the state space model (SSM)
#  expectations.
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Revision History
#   Wed Nov 14 13:34:01 Central Standard Time 2012 -- Michael Hunter created file
#   Sat Nov 17 16:14:56 Central Standard Time 2012 -- Michael Hunter change names to ExpectationStateSpace instead of ExpectationSSM
#   Wed Sep 18 23:09:34 Central Daylight Time 2013 -- Michael Hunter changed name of arguments to x0 and P0
#   

#--------------------------------------------------------------------
# **DONE**
setClass(Class = "MxExpectationStateSpace",
	representation = representation(
		A = "MxCharOrNumber",
		B = "MxCharOrNumber",
		C = "MxCharOrNumber",
		D = "MxCharOrNumber",
		Q = "MxCharOrNumber",
		R = "MxCharOrNumber",
		x0 = "MxCharOrNumber",
		P0 = "MxCharOrNumber",
		u = "MxCharOrNumber",
		dims = "character",
		definitionVars = "list",
		t = "MxCharOrNumber",
		scores = "logical",
		AIsZero = "logical",
		xPredicted = "matrix",
		yPredicted = "matrix",
		PPredicted = "matrix",
		SPredicted = "matrix",
		xUpdated = "matrix",
		PUpdated = "matrix",
		xSmoothed = "matrix",
		PSmoothed = "matrix"),
	contains = "MxBaseExpectation")


#--------------------------------------------------------------------
# **DONE**
setMethod("initialize", "MxExpectationStateSpace",
	function(.Object, A, B, C, D, Q, R, x0, P0, u, dims, t, scores,
		data = as.integer(NA), name = 'expectation') {
		.Object@name <- name
		.Object@A <- A
		.Object@B <- B
		.Object@C <- C
		.Object@D <- D
		.Object@Q <- Q
		.Object@R <- R
		.Object@x0 <- x0
		.Object@P0 <- P0
		.Object@u <- u
		.Object@data <- data
		.Object@dims <- dims
		.Object@t <- t
		.Object@scores <- scores
		.Object@AIsZero <- FALSE
		.Object@definitionVars <- list()
		return(.Object)
	}
)


#--------------------------------------------------------------------
# TODO: Ask Spiegel and Brick what this is really supposed to do.
setMethod("genericExpConvertEntities", "MxExpectationStateSpace",
	function(.Object, flatModel, namespace, labelsData) {
		cache <- new.env(parent = emptyenv())
		modelname <- getModelName(.Object)
		if(is.na(.Object@data)) {
			msg <- paste("The SSM expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}
		
		tuple <- evaluateMxObject(.Object@A, flatModel, labelsData, cache)
		Amatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@B, flatModel, labelsData, cache)
		Bmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@C, flatModel, labelsData, cache)
		Cmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@D, flatModel, labelsData, cache)
		Dmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@Q, flatModel, labelsData, cache)
		Qmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@R, flatModel, labelsData, cache)
		Rmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@x0, flatModel, labelsData, cache)
		Xmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@P0, flatModel, labelsData, cache)
		Pmatrix <- tuple[[1]]
		cache <- tuple[[2]]
		tuple <- evaluateMxObject(.Object@u, flatModel, labelsData, cache)
		Umatrix <- tuple[[1]]
		cache <- tuple[[2]]
		if(!is.na(.Object@t)){
			tuple <- evaluateMxObject(.Object@t, flatModel, labelsData, cache)
			Tmatrix <- tuple[[1]]
			cache <- tuple[[2]]
		} else {Tmatrix <- NULL}
		# Conformability checking
		ldim <- ncol(Cmatrix)
		mdim <- nrow(Cmatrix)
		checkSSMConformable(Amatrix, ldim, ldim, 'A', omxQuotes(modelname))
		checkSSMConformable(Cmatrix, mdim, ldim, 'C', omxQuotes(modelname))
		checkSSMConformable(Qmatrix, ldim, ldim, 'Q', omxQuotes(modelname))
		checkSSMConformable(Rmatrix, mdim, mdim, 'R', omxQuotes(modelname))
		checkSSMConformable(Xmatrix, ldim, 1, 'x0', omxQuotes(modelname))
		checkSSMConformable(Pmatrix, ldim, ldim, 'P0', omxQuotes(modelname))
		#
		# If any of B, D, u are not missing, then
		#  1.  All of B, D, u must not be missing
		#  2.  colnames of B and D must match rownames of u?
		# 3.   Conformability checks
		#  ldim <- ncol(Cmatrix)
		#  mdim <- nrow(Cmatrix)
		#  udim <- ncol(Bmatrix)
		#  A: ldim x ldim
		#  B: ldim x udim
		#  C: mdim x ldim
		#  D: mdim x udim
		#  Q: ldim x ldim
		#  R: mdim x mdim
		#  x: ldim x 1
		#  P: ldim x ldim
		#  u: udim x 1
		if(!is.null(Bmatrix) || !is.null(Dmatrix) || !is.null(Umatrix)){
			udim <- ncol(Bmatrix)
			checkSSMConformable(Bmatrix, ldim, udim, 'B', omxQuotes(modelname))
			checkSSMConformable(Dmatrix, mdim, udim, 'D', omxQuotes(modelname))
			checkSSMConformable(Umatrix, udim, 1, 'u', omxQuotes(modelname))
		}
		flatModel <- updateSSMdimnames(.Object, flatModel)
		
		return(flatModel)
	}
)

updateSSMdimnames <- function(flatExpectation, flatJob){
	cMatrixName <- flatExpectation@C
	cMatrix <- flatJob[[cMatrixName]]
	if (is.null(cMatrix)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown C matrix name", 
			omxQuotes(simplifyName(cMatrixName, modelname)),
			"detected in the State Space expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatExpectation@dims
	if (!is.null(dimnames(cMatrix)) && !single.na(dims) && 
		!identical(dimnames(cMatrix)[[1]], dims)) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The C matrix associated",
			"with the State Space expectation function in model", 
			omxQuotes(modelname), "contains dimnames and",
			"the expectation function has specified different dimnames")
		stop(msg, call.=FALSE)		
	}
	if (is.null(dimnames(cMatrix)) && !single.na(dims)) {
		dimnames(flatJob[[cMatrixName]]) <- list(dims, c())
	}
	
	return(flatJob)
}


#--------------------------------------------------------------------
# **DONE**
setMethod("qualifyNames", signature("MxExpectationStateSpace"), 
	function(.Object, modelname, namespace) {
		.Object@name <- imxIdentifier(modelname, .Object@name)
		.Object@A <- imxConvertIdentifier(.Object@A, modelname, namespace)
		.Object@B <- imxConvertIdentifier(.Object@B, modelname, namespace)
		.Object@C <- imxConvertIdentifier(.Object@C, modelname, namespace)
		.Object@D <- imxConvertIdentifier(.Object@D, modelname, namespace)
		.Object@Q <- imxConvertIdentifier(.Object@Q, modelname, namespace)
		.Object@R <- imxConvertIdentifier(.Object@R, modelname, namespace)
		.Object@x0 <- imxConvertIdentifier(.Object@x0, modelname, namespace)
		.Object@P0 <- imxConvertIdentifier(.Object@P0, modelname, namespace)
		.Object@u <- imxConvertIdentifier(.Object@u, modelname, namespace)
		.Object@data <- imxConvertIdentifier(.Object@data, modelname, namespace)
		.Object@t <- imxConvertIdentifier(.Object@t, modelname, namespace)
		#.Object@scores <- imxConvertIdentifier(.Object@scores, modelname, namespace)
		return(.Object)
	}
)




#--------------------------------------------------------------------

checkSSMNotMissing <- function(matrixobj, matrixname, modelname){
	if(is.null(matrixobj)){
		msg <- paste("The", matrixname, "matrix of the state space expectation in model",
			modelname, "is missing.")
		stop(msg, call. = FALSE)
	}
}

checkSSMConformable <- function(mat, rows, cols, matname, modname){
	if( nrow(mat) != rows || ncol(mat) != cols ){
		msg <- paste("The ", matname, " matrix is not the correct size",
			" in the state space expectation of model ", modname,
			".  It is ", nrow(mat), " by ", ncol(mat), " and should be ",
			rows, " by ", cols, ".", sep="")
		stop(msg, call. = FALSE)
	}
}

# TODO: Allow subsets of the matrices to be specified
#  by filling in default matrices.
setMethod("genericExpFunConvert", signature("MxExpectationStateSpace"), 
	function(.Object, flatModel, model, labelsData, dependencies) {
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]	
		name <- .Object@name
		aMatrix <- .Object@A
		bMatrix <- .Object@B
		cMatrix <- .Object@C
		dMatrix <- .Object@D
		qMatrix <- .Object@Q
		rMatrix <- .Object@R
		xMatrix <- .Object@x0
		pMatrix <- .Object@P0
		uMatrix <- .Object@u
		tMatrix <- .Object@t
		data <- .Object@data
		if(is.na(data)) {
			msg <- paste("The state space expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call. = FALSE)
		}
		mxDataObject <- flatModel@datasets[[.Object@data]]
		checkNumericData(mxDataObject)
		.Object@A <- imxLocateIndex(flatModel, aMatrix, name)
		.Object@B <- imxLocateIndex(flatModel, bMatrix, name)
		.Object@C <- imxLocateIndex(flatModel, cMatrix, name)
		.Object@D <- imxLocateIndex(flatModel, dMatrix, name)
		.Object@Q <- imxLocateIndex(flatModel, qMatrix, name)
		.Object@R <- imxLocateIndex(flatModel, rMatrix, name)
		.Object@x0 <- imxLocateIndex(flatModel, xMatrix, name)
		.Object@P0 <- imxLocateIndex(flatModel, pMatrix, name)
		.Object@u <- imxLocateIndex(flatModel, uMatrix, name)
		.Object@t <- imxLocateIndex(flatModel, tMatrix, name)
		.Object@data <- as.integer(imxLocateIndex(flatModel, data, name))
		#
		# Check the data has row and column names as appropriate
		verifyObservedNames(mxDataObject@observed, mxDataObject@means, mxDataObject@type, flatModel, modelname, "State Space")
		#
		# Change *Matrix from the string name of the matrix to the object
		aMatrix <- flatModel[[aMatrix]]
		bMatrix <- flatModel[[bMatrix]]
		cMatrix <- flatModel[[cMatrix]]
		dMatrix <- flatModel[[dMatrix]]
		qMatrix <- flatModel[[qMatrix]]
		rMatrix <- flatModel[[rMatrix]]
		xMatrix <- flatModel[[xMatrix]]
		pMatrix <- flatModel[[pMatrix]]
		uMatrix <- flatModel[[uMatrix]]
		tMatrix <- flatModel[[tMatrix]]
		#
		# Check for missing matrices or combinations of matrices
		checkSSMNotMissing(cMatrix, 'C', omxQuotes(modelname))
		# Check C for dimnames
		if (is.null(dimnames(cMatrix))) {
			msg <- paste("The C matrix of model",
				omxQuotes(modelname), "does not contain dimnames")
			stop(msg, call. = FALSE)
		}
		if (is.null(dimnames(cMatrix)[[2]])) {
			msg <- paste("The C matrix of model",
				omxQuotes(modelname), "does not contain colnames")
			stop(msg, call. = FALSE)
		}
		checkSSMNotMissing(aMatrix, 'A', omxQuotes(modelname))
		checkSSMNotMissing(qMatrix, 'Q', omxQuotes(modelname))
		checkSSMNotMissing(rMatrix, 'R', omxQuotes(modelname))
		checkSSMNotMissing(xMatrix, 'x0', omxQuotes(modelname))
		checkSSMNotMissing(pMatrix, 'P0', omxQuotes(modelname))
		if(!is.null(bMatrix) || !is.null(dMatrix) || !is.null(uMatrix)){
			checkSSMNotMissing(bMatrix, 'B', omxQuotes(modelname))
			checkSSMNotMissing(dMatrix, 'D', omxQuotes(modelname))
			checkSSMNotMissing(uMatrix, 'u', omxQuotes(modelname))
		}
		#
		# If A is a matrix (not algebra), has all zero values, all NA labels, no free parameters
		if("MxMatrix" %in% is(aMatrix) && all(aMatrix$values == 0) && all(is.na(aMatrix$labels)) && all(aMatrix$free == FALSE)){
			.Object@AIsZero <- TRUE
			# This is used in the backend for continuous time models to allow for constant slope models
		}
		#
		# Do data processing
		translatedNames <- c(dimnames(cMatrix)[[1]])
		if (mxDataObject@type == 'raw') {
			.Object@dataColumnNames <- translatedNames
			.Object@dataColumns <- generateDataColumns(flatModel, translatedNames, data)
			if (length(mxDataObject@observed) == 0) {
				.Object@data <- as.integer(NA)
			}
			if (single.na(.Object@dims)) {
				.Object@dims <- translatedNames
			}
		}
		return(.Object)
	}
)


#--------------------------------------------------------------------
# **DONE**
setMethod("genericExpDependencies", signature("MxExpectationStateSpace"),
	function(.Object, dependencies) {
		sources <- c(.Object@A, .Object@B, .Object@C, .Object@D, .Object@Q, .Object@R,
			.Object@x0, .Object@P0, .Object@u, .Object@t)
		sources <- sources[!is.na(sources)]
		dependencies <- imxAddDependency(sources, .Object@name, dependencies)
		return(dependencies)
	}
)


#--------------------------------------------------------------------
# **DONE**
setMethod("genericExpRename", signature("MxExpectationStateSpace"),
	function(.Object, oldname, newname) {
		.Object@A <- renameReference(.Object@A, oldname, newname)
		.Object@B <- renameReference(.Object@B, oldname, newname)
		.Object@C <- renameReference(.Object@C, oldname, newname)
		.Object@D <- renameReference(.Object@D, oldname, newname)
		.Object@Q <- renameReference(.Object@Q, oldname, newname)
		.Object@R <- renameReference(.Object@R, oldname, newname)
		.Object@x0 <- renameReference(.Object@x0, oldname, newname)
		.Object@P0 <- renameReference(.Object@P0, oldname, newname)
		.Object@u <- renameReference(.Object@u, oldname, newname)
		.Object@data <- renameReference(.Object@data, oldname, newname)
		.Object@t <- renameReference(.Object@t, oldname, newname)
		return(.Object)
	}
)


getAllModelNames <- function(model){
	ret <- getModelName(model)
	if(length(model@submodels) > 0){
		ret <- c(ret, unlist(lapply(model@submodels, getAllModelNames)))
	}
	return(ret)
}

#--------------------------------------------------------------------

setMethod("genericGetExpected", signature("MxExpectationStateSpace"),
		function(.Object, model, what, defvar.row=1, subname=model@name) {
			ret <- list()
			if(length(defvar.row) > 1){
				stop("'defvar.row' must be (1) a single integer, (2) 'all', or (3) Inf")
			}
			wantMean <- any(c('mean', 'means') %in% what)
			wantCov <- any(c('covariance','covariances') %in% what)
			if(defvar.row == Inf){
				if(!single.na(.Object@t)){
					stop("Found continuous time model.\nAsymptotic expectations are not yet implemented for continuous time models.")
				}
				if(imxHasDefinitionVariable(model)){
					stop("Found definition variables.\nAsymptotic expectations are not valid for models with definition variables.")
				}
				#stop("Frau Bl\u00FCcher! This is not yet implemented.")
				Aname <- paste(subname, .Object@A, sep=".")
				Qname <- paste(subname, .Object@Q, sep=".")
				Cname <- paste(subname, .Object@C, sep=".")
				Rname <- paste(subname, .Object@R, sep=".")
				A <- mxEvalByName(Aname, model, compute=TRUE)
				Q <- mxEvalByName(Qname, model, compute=TRUE)
				C <- mxEvalByName(Cname, model, compute=TRUE)
				R <- mxEvalByName(Rname, model, compute=TRUE)
				I <- diag(1, nrow=nrow(A)*nrow(A))
				ImA <- try(solve(I - A %x% A))
				if(class(ImA) %in% "try-error"){
					stop("Could not invert I - A %x% A\nAsymptotic expectations are not valid in this case.")
				}
				Pinf <- ImA %*% matrix(c(Q), ncol=1)
				Pinf <- matrix(Pinf, nrow=nrow(A), ncol=nrow(A))
				Sinf <- C %*% Pinf %*% t(C) + R
				if (wantCov) {
					ret[['covariance']] <- Sinf
				}
				if (wantMean) {
					ret[['means']] <- 'Not implemented'
				}
			} else {
				ks <- mxKalmanScores(model, frontend=FALSE)
				if(defvar.row =="all"){
					defvar.row <- 1:(nrow(ks$xPredicted)-1)
				} else {
					defvar.row <- defvar.row + 1
				}
				if (wantCov) {
					ret[['covariance']] <- ks$SPredicted[ , , defvar.row, drop=FALSE]
				}
				if (wantMean) {
					ret[['means']] <- ks$yPredicted[defvar.row, , drop=FALSE]
				}
			}
			ret
})


#--------------------------------------------------------------------
checkSSMargument <- function(x, xname) {
	if (!(single.na(x) || typeof(x) == "character")) {
		msg <- paste("argument ", xname, " is not a string ",
			"(the name of the '", xname, "' matrix)", sep="")
		stop(msg)
	}
	if (is.na(x)) x <- as.integer(NA)
	return(x)
}

#--------------------------------------------------------------------
# **DONE**
mxExpectationStateSpace <- function(A, B, C, D, Q, R, x0, P0, u, dimnames = NA, thresholds = deprecated(), threshnames = deprecated(), ..., t=NA, scores=FALSE){
	A <- checkSSMargument(A, "A")
	B <- checkSSMargument(B, "B")
	C <- checkSSMargument(C, "C")
	D <- checkSSMargument(D, "D")
	Q <- checkSSMargument(Q, "Q")
	R <- checkSSMargument(R, "R")
	x0 <- checkSSMargument(x0, "x0")
	P0 <- checkSSMargument(P0, "P0")
	u <- checkSSMargument(u, "u")
	t <- checkSSMargument(t, "t")
	if (length(scores) > 1 || typeof(scores) != "logical") {
		stop("'scores' argument is not a logical value")
	}
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("Dimnames argument is not a character vector")
	}
	if (length(dimnames) == 0) {
		stop("Dimnames argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for dimnames vector")
	}
	return(new("MxExpectationStateSpace", A, B, C, D, Q, R, x0, P0, u, dimnames, t=t, scores=scores))
}

mxExpectationStateSpaceContinuousTime <- function(A, B, C, D, Q, R, x0, P0, u, t=NA, dimnames = NA, thresholds = deprecated(), threshnames = deprecated(), ..., scores=FALSE){
	mxExpectationStateSpace(t=t, scores=scores, A, B, C, D, Q, R, x0, P0, u, dimnames)
}

mxExpectationSSCT <- mxExpectationStateSpaceContinuousTime


#--------------------------------------------------------------------
# TODO: Add expected mean and cov printouts
displayMxExpectationStateSpace <- function(expectation) {
	cat("MxExpectationStateSpace", omxQuotes(expectation@name), '\n')
	cat("$A :", omxQuotes(expectation@A), '\n')
	cat("$B :", omxQuotes(expectation@B), '\n')
	cat("$C :", omxQuotes(expectation@C), '\n')
	cat("$D :", omxQuotes(expectation@D), '\n')
	cat("$Q :", omxQuotes(expectation@Q), '\n')
	cat("$R :", omxQuotes(expectation@R), '\n')
	cat("$x0 :", omxQuotes(expectation@x0), '\n')
	cat("$P0 :", omxQuotes(expectation@P0), '\n')
	cat("$u :", omxQuotes(expectation@u), '\n')
	if (single.na(expectation@dims)) {
		cat("$dims : NA \n")
	} else {
		cat("$dims :", omxQuotes(expectation@dims), '\n')
	}		
	invisible(expectation)
}

setMethod("print", "MxExpectationStateSpace", function(x,...) { 
	displayMxExpectationStateSpace(x) 
})

setMethod("show", "MxExpectationStateSpace", function(object) { 
	displayMxExpectationStateSpace(object) 
})


#--------------------------------------------------------------------
KalmanFilter <- function(A, B, C, D, Q, R, x, y, u, P, ct=FALSE, dt=0){
	if(ct){
		I <- diag(1, nrow=nrow(A))
		out <- kalmanDiscretize(A, B, Q, dt, I)
		A <- out$Ad
		B <- out$Bd
		Q <- out$Qd
	}
	
	x <- A %*% x + B %*% u
	P <- A %*% P %*% t(A) + Q
	x.pred <- x
	P.pred <- P
	
	r <- y - (C %*% x + D %*% u)
	notMiss <- !is.na(r)
	r[!notMiss] <- 0
	if(length(r)==sum(!notMiss)){#all missing row
		m2ll <- log(det(C %*% P %*% t(C) + R))
		return(list(x.pred=x.pred, P.pred=P.pred, x.upda=x.pred, P.upda=P.pred, m2ll=m2ll, L=exp(m2ll/-2), A=A ))
	} else {
		Cf <- C[notMiss, , drop=FALSE]
		Rf <- R[notMiss, notMiss, drop=FALSE]
		S <- Cf %*% P %*% t(Cf) + Rf
		Sinv <- solve(S)
		rf <- matrix(r[notMiss], ncol=1)
		K <- P %*% t(Cf) %*% Sinv
		x <- x + K %*% rf
		P <- P - K %*% Cf %*% P
		x.upda <- x
		P.upda <- P
		
		const <- length(rf)*log(2*pi)
		m2ll <- log(det(S)) + t(rf) %*% Sinv %*% rf + const
		
		return(list(x.pred=x.pred, P.pred=P.pred, x.upda=x.upda, P.upda=P.upda, m2ll=m2ll, L=exp(m2ll/-2), A=A ))
	}
}

kalmanBackendScoreHelper <- function(model, data=NA){
	if(!single.na(data)){
		model@data <- mxData(data, type='raw')
	}
	param <- coef(model)
	if(length(param) > 0){
		model <- omxSetParameters(model, labels=names(param), free=FALSE)
	}
	e <- model@expectation
	model <- mxModel(model=model, name='KalmanScoring',
		mxExpectationStateSpace(
			A=e@A, B=e@B, C=e@C, D=e@D,
			Q=e@Q, R=e@R, x0=e@x0, P0=e@P0, u=e@u,
			t=e@t, scores=TRUE),
		mxFitFunctionML(rowDiagnostics=TRUE))
	run <- mxRun(model, silent=TRUE)
	re <- run@expectation
	fitfunction <- NULL
	L <- c(1, attr(mxEval(fitfunction, run), 'likelihoods'))
	sp <- vech2array(re@SPredicted)
	sp[,,1] <- NA
	yp <- re@yPredicted
	yp[1,] <- NA
	
	return(list(xPredicted=re@xPredicted, PPredicted=vech2array(re@PPredicted), xUpdated=re@xUpdated, PUpdated=vech2array(re@PUpdated), xSmoothed=re@xSmoothed, PSmoothed=vech2array(re@PSmoothed), m2ll=-2*log(L), L=L, yPredicted=yp, SPredicted=sp))
}

vech2array <- function(x){
	y <- apply(x, 1, vech2full)
	xdim <- ifelse(is.matrix(y), sqrt(nrow(y)), 1)
	tdim <- ifelse(is.matrix(y), ncol(y), length(y))
	array(y, c(xdim, xdim, tdim))
}

kalmanFrontendScoreHelper <- function(model, data=NA){
	message("Computing Kalman scores in frontend R.  This may take a few seconds.")
	if(single.na(data)) {
		#TODO check that data are raw
		data <- model@data@observed
	}
	x0 <- mxEvalByName(model@expectation@x0, model, compute=TRUE)
	P0 <- mxEvalByName(model@expectation@P0, model, compute=TRUE)
	A <- mxEvalByName(model@expectation@A, model, compute=TRUE)
	B <- mxEvalByName(model@expectation@B, model, compute=TRUE)
	C <- mxEvalByName(model@expectation@C, model, compute=TRUE)
	D <- mxEvalByName(model@expectation@D, model, compute=TRUE)
	Q <- mxEvalByName(model@expectation@Q, model, compute=TRUE)
	R <- mxEvalByName(model@expectation@R, model, compute=TRUE)
	u <- mxEvalByName(model@expectation@u, model, compute=TRUE)
	
	
	hasDefVars <- imxHasDefinitionVariable(model)
	continuousTime <- !single.na(model@expectation@t)
	
	X.pred <- matrix(0, nrow=nrow(data)+1, ncol=nrow(x0))
	X.upda <- matrix(0, nrow=nrow(data)+1, ncol=nrow(x0))
	X.pred[1,] <- x0
	X.upda[1,] <- x0
	P.pred <- array(0, dim=c(nrow(x0), nrow(x0), nrow(data)+1))
	P.upda <- array(0, dim=c(nrow(x0), nrow(x0), nrow(data)+1))
	A.disc <- array(0, dim=c(nrow(x0), nrow(x0), nrow(data)))
	P.pred[,,1] <- P0
	P.upda[,,1] <- P0
	m2ll <- numeric(nrow(data)+1)
	m2ll[1] <- 0
	L <- numeric(nrow(data)+1)
	L[1] <- 1
	oldT <- 0
	newT <- 0
	for(i in 1:nrow(data)){
		if(hasDefVars){
			A <- mxEvalByName(model@expectation@A, model, compute=TRUE, defvar.row=i)
			B <- mxEvalByName(model@expectation@B, model, compute=TRUE, defvar.row=i)
			C <- mxEvalByName(model@expectation@C, model, compute=TRUE, defvar.row=i)
			D <- mxEvalByName(model@expectation@D, model, compute=TRUE, defvar.row=i)
			Q <- mxEvalByName(model@expectation@Q, model, compute=TRUE, defvar.row=i)
			R <- mxEvalByName(model@expectation@R, model, compute=TRUE, defvar.row=i)
			u <- mxEvalByName(model@expectation@u, model, compute=TRUE, defvar.row=i)
			if(continuousTime){
				newT <- mxEvalByName(model@expectation@t, model, compute=TRUE, defvar.row=i)
			}
		}
		
		deltaT <- newT - oldT
		oldT <- newT
		
		res <- KalmanFilter(A=A, B=B, C=C, D=D, Q=Q, R=R, x=matrix(X.upda[i,]), y=matrix(unlist(data[i,rownames(C)])), u=u, P=P.upda[,,i], continuousTime, deltaT)
		
		X.pred[i+1,] <- res$x.pred
		X.upda[i+1,] <- res$x.upda
		P.pred[,,i+1] <- res$P.pred
		P.upda[,,i+1] <- res$P.upda
		m2ll[i+1] <- res$m2ll
		L[i+1] <- res$L
		A.disc[,,i] <- res$A
	}
	X.smoo <- matrix(0, nrow=nrow(data)+1, ncol=nrow(x0))
	X.smoo[nrow(data)+1,] <- X.upda[nrow(data)+1, ]
	P.smoo <- array(0, dim=c(nrow(x0), nrow(x0), nrow(data)+1))
	P.smoo[,,nrow(data)+1] <- P.upda[,,nrow(data)+1]
	for(i in nrow(data):1){
		A <- A.disc[,,i]
		SGain <- P.upda[,,i] %*% A %*% solve(P.pred[,,i+1])
		X.smoo[i,] <- matrix(X.upda[i,]) + SGain %*% matrix(X.smoo[i+1,] - X.pred[i+1,])
		P.smoo[,,i] <- P.upda[,,i] + SGain %*% (P.smoo[,,i+1] - P.pred[,,i+1]) %*% t(SGain)
	}
	
	return(list(xPredicted=X.pred, PPredicted=P.pred, xUpdated=X.upda, PUpdated=P.upda, xSmoothed=X.smoo, PSmoothed=P.smoo, m2ll=m2ll, L=L))
}

mxKalmanScores <- function(model, data=NA, frontend=TRUE){
  warnModelCreatedByOldVersion(model)
	if(!frontend){
		scores <- kalmanBackendScoreHelper(model, data)
	} else {
		scores <- kalmanFrontendScoreHelper(model, data)
	}
	return(scores)
}


# From StateSpaceOsc.R
#smod2 <- mxModel(name='with scores', model=smod,
#	mxExpectationStateSpace(A='A', B='B', C='C', D='D', Q='Q', R='R', x0='x', P0='P', u='u', scores=TRUE)
#)
#smod3 <- omxSetParameters(smod2, labels=names(omxGetParameters(smod2)), free=FALSE)
#srun3 <- mxRun(smod3)
#a <- Sys.time(); res <- mxKalmanScores(srun); b <- Sys.time()
#as.numeric(b-a)/as.numeric(summary(srun3)$wallTime)
## [1] 147.7429
# Without the backward pass yet (i.e. no smoother)
#  backend is 147 times as fast as frontent.
#require(rbenchmark)
#benchmark(mxKalmanScores(srun), mxRun(smod3), replications=20)
## relative time: 186x speedup
#benchmark(mxRun(smod), mxRun(smod2), replications=20)
## .2% increase in relative time
## .0085 extra seconds per model

#res <- mxKalmanScores(srun)


#--------------------------------------------------------------------
setMethod("genericGenerateData", signature("MxExpectationStateSpace"),
	function(.Object, model, nrows, subname, empirical, returnModel, use.miss,
		   .backend, nrowsProportion) {
  origData <- findDataForSubmodel(model, subname)
  origRows <- if (!is.null(origData)) { nrowMxData(origData) } else { NULL }
  nrows <- calcNumRows(nrows, nrowsProportion, origRows, subname)

		A <- mxEvalByName(model[[subname]]@expectation@A, model, compute=TRUE)
		B <- mxEvalByName(model[[subname]]@expectation@B, model, compute=TRUE)
		C <- mxEvalByName(model[[subname]]@expectation@C, model, compute=TRUE)
		D <- mxEvalByName(model[[subname]]@expectation@D, model, compute=TRUE)
		Q <- mxEvalByName(model[[subname]]@expectation@Q, model, compute=TRUE)
		R <- mxEvalByName(model[[subname]]@expectation@R, model, compute=TRUE)
		u <- mxEvalByName(model[[subname]]@expectation@u, model, compute=TRUE)
		
		hasDefVars <- imxHasDefinitionVariable(model)
		continuousTime <- !single.na(model[[subname]]@expectation@t)
		
		x0 <- mxEvalByName(model[[subname]]@expectation@x0, model, compute=TRUE)
		P0 <- mxEvalByName(model[[subname]]@expectation@P0, model, compute=TRUE)
		
		tdim <- nrows
		ydim <- nrow(C)
		xdim <- nrow(A)
		tx <- matrix(0, xdim, tdim+1)
		ty <- matrix(0, ydim, tdim)
		I <- diag(1, nrow=nrow(A))
		Z <- diag(0, nrow=nrow(A))
		
		tx[,1] <- x0
		oldT <- 0
		for(i in 2:(tdim+1)){
			if(hasDefVars){
				A <- mxEvalByName(model[[subname]]@expectation@A, model, compute=TRUE, defvar.row=i-1)
				B <- mxEvalByName(model[[subname]]@expectation@B, model, compute=TRUE, defvar.row=i-1)
				C <- mxEvalByName(model[[subname]]@expectation@C, model, compute=TRUE, defvar.row=i-1)
				D <- mxEvalByName(model[[subname]]@expectation@D, model, compute=TRUE, defvar.row=i-1)
				Q <- mxEvalByName(model[[subname]]@expectation@Q, model, compute=TRUE, defvar.row=i-1)
				R <- mxEvalByName(model[[subname]]@expectation@R, model, compute=TRUE, defvar.row=i-1)
				u <- mxEvalByName(model[[subname]]@expectation@u, model, compute=TRUE, defvar.row=i-1)
				newT <- mxEvalByName(model[[subname]]@expectation@t, model, compute=TRUE, defvar.row=i-1)
			}
			if(continuousTime){
				#browser()
				deltaT <- c(newT - oldT)
				oldT <- newT
				kd <- kalmanDiscretize(A, B, Q, deltaT, I)
				Ad <- kd$Ad
				Bd <- kd$Bd
				Qd <- kd$Qd
			} else {
				Ad <- A
				Bd <- B
				Qd <- Q
			}
			xp <- Ad %*% tx[,i-1] + Bd %*% u
			tx[,i] <- xp + t(.rmvnorm(1, rep(0, xdim), Qd, empirical))
			ty[,i-1] <- C %*% tx[,i-1] + D %*% u + t(.rmvnorm(1, rep(0, ydim), R, empirical))
		}
		ret <- t(ty)
		colnames(ret) <- dimnames(C)[[1]]
		if (returnModel) {
		  mxModel(model[[subname]], mxData(as.data.frame(ret), "raw"))
		} else {
		  as.data.frame(ret)
		}
	}
)

kalmanDiscretize <- function(A, B, Q, deltaT, I){
	if(length(deltaT) > 1){
		stop('Found bad time argument for continuous time state space model.')
	}
	deltaT <- c(deltaT)
	
	BLOCK <- matrix(0, nrow=2*nrow(A), ncol=2*ncol(A))
	
	# First Block expm for A integral, and expm(A*deltaT)
	BLOCK[1:(2*nrow(A)), 1:ncol(A)] <- 0
	BLOCK[1:nrow(A), (nrow(A)+1):(2*nrow(A))] <- I
	BLOCK[(nrow(A)+1):(2*nrow(A)), (nrow(A)+1):(2*nrow(A))] <- A
	BLOCK <- OpenMx::expm(BLOCK*deltaT)
	expA <- BLOCK[(nrow(A)+1):(2*nrow(A)), (nrow(A)+1):(2*nrow(A))]
	intA <- BLOCK[1:nrow(A), (nrow(A)+1):(2*nrow(A))]
	
	# Second Block expm for discretized Q
	BLOCK[1:(nrow(A)), 1:ncol(A)] <- -t(A)
	BLOCK[(nrow(A)+1):(2*nrow(A)), 1:ncol(A)] <- 0
	BLOCK[1:nrow(A), (nrow(A)+1):(2*nrow(A))] <- Q
	BLOCK[(nrow(A)+1):(2*nrow(A)), (nrow(A)+1):(2*nrow(A))] <- A
	BLOCK <- OpenMx::expm(BLOCK*deltaT)
	
	Bd <- intA %*% B
	Qd <- t(expA) %*% BLOCK[1:nrow(A), (nrow(A)+1):(2*nrow(A))]
	return(list(Ad=expA, Bd=Bd, Qd=Qd))
}

