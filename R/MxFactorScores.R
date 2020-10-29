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
# Author: Michael D. Hunter
# Date: 2015-06-14
# Filename: MxFactorScores.R
# Purpose: Write a helper function for computing various type of factor scores
#------------------------------------------------------------------------------


requireMinManifests <- function(row) {
	stop(paste("mxFactorScores: row", row, "has missing data.",
		   "Hence, you must specify minManifests"),
	     call. = FALSE)
}

mxFactorScores <- function(model, type=c('ML', 'WeightedML', 'Regression'), minManifests=as.integer(NA))
{
  warnModelCreatedByOldVersion(model)
	if(length(unlist(strsplit(model@name, split=' ', fixed=TRUE))) > 1){
		message(paste('The model called', omxQuotes(model@name), 'has spaces in the model name.  I cannot handle models with spaces in the model name, so I removed them before getting factor scores.'))
		model <- mxRename(model, paste(unlist(strsplit(model@name, split=' ', fixed=TRUE)), collapse=''))
	}
	# Handling of multigroup models
	if(is.null(model$expectation) && (class(model$fitfunction) %in% "MxFitFunctionMultigroup") ){
		submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
		ret <- list()
		for(amod in submNames){
			ret[[amod]] <- mxFactorScores(model[[amod]], type, minManifests)
		}
		return(ret)
	}
	if(model$data$type!='raw'){
		stop("The 'model' argument must have raw (not summary) data.")
	}
	classExpect <- class(model$expectation)
	if(!(classExpect %in% "MxExpectationLISREL") && !(classExpect %in% "MxExpectationRAM")){
		stop('Factor scores are only implemented for LISREL and RAM expectations.')
	}
	if((classExpect %in% "MxExpectationLISREL") && !single.na(model$expectation$LY)){
		stop('Factor scores for LISREL are only implemented for the exogenous-only model, but an LY matrix was detected.  Try restructuring your model as LISREL exogenous-only, or as RAM.')
	}
	if(classExpect %in% "MxExpectationLISREL"){
		lx <- mxEvalByName(model$expectation$LX, model, compute=TRUE)
		nksix <- dim(lx)
		nksi <- nksix[2]
		nx <- nksix[1]
		factorNames <- dimnames(lx)[[2]]
		factorScoreHelperFUN <- lisrelFactorScoreHelper
	} else if(classExpect %in% "MxExpectationRAM"){
		fm <- mxEvalByName(model$expectation$F, model, compute=TRUE)
		nksix <- dim(fm)
		nksi <- nksix[2] - nksix[1]
		nx <- nksix[1]
		factorNames <- dimnames(fm)[[2]][!createOppositeF(fm)$is.manifest]
		factorScoreHelperFUN <- ramFactorScoreHelper
	}
	nrows <- nrow(model$data$observed)
	res <- array(as.numeric(NA), c(nrows, nksi, 2))
	if(any(type %in% c('ML', 'WeightedML'))){
		model <- omxSetParameters(model, labels=names(omxGetParameters(model)), free=FALSE)
		work <- factorScoreHelperFUN(model)
		dataModelName <- work$name
		if(type[1]=='WeightedML'){
			wup <- mxModel(model="Container", work,
				mxAlgebraFromString(paste(work@name, ".weight + ", work@name, ".fitfunction", sep=""), name="wtf"),
				mxFitFunctionAlgebra("wtf")
			)
			work <- wup
		}
		work@data <- NULL
		fullData <- as.data.frame(model$data$observed)
		plan1 <- list(
			GD=mxComputeGradientDescent(nudgeZeroStarts=FALSE))
		wantSE <- tolower(mxOption(model=model, key="Standard Errors")) == "yes"
		if (wantSE) {
			plan1 <- c(plan1,
				ND=mxComputeNumericDeriv(),
				SE=mxComputeStandardError())
		}
		plan <- list(
			SOS=mxComputeSetOriginalStarts(),
			LD=mxComputeLoadData(dataModelName, names(fullData),
				method="data.frame", byrow=FALSE,
				observed=fullData),
			TC=mxComputeTryCatch(mxComputeSequence(plan1)),
			CP=mxComputeCheckpoint(toReturn=TRUE, standardErrors = wantSE))
		plan <- mxComputeLoop(plan, maxIter = nrow(fullData))
		rawData <- fullData[1,,drop=FALSE]
		if(type[1]=='ML'){
			fit <- mxModel(model=work, mxData(rawData, 'raw'))
		} else if(type[1]=='WeightedML'){
			work@submodels[[1]]@data <- mxData(rawData, 'raw')
			fit <- mxModel(model=work)
		}
		fit <- mxRun(mxModel(fit, plan))
		got <- fit$compute$steps$CP$log
		res[,,1] <- as.matrix(got[,names(coef(fit)),drop=FALSE])
		if (wantSE) {
			res[,,2] <- as.matrix(got[,paste0(names(coef(fit)), 'SE'),drop=FALSE])
		} else {
			msg <- paste0("factor-score standard errors not available from MxModel '",
				model$name,"' because calculating SEs is turned off for that ",
				"model (possibly due to one or more MxConstraints)")
			warning(msg, sep="")
		}
		res <- clearExcessivelyMissingRows(fullData, minManifests, res)
	} else if(tolower(type)=='regression'){
		if(!single.na(model$expectation$thresholds)){
			stop('Regression factor scores cannot be computed when there are thresholds (ordinal data).')
		}
		if(!(classExpect %in% "MxExpectationLISREL")){
			#stop('Regression factor scores are only possible for LISREL expectations.')
			res <- RAMrfs(model, res, minManifests)
		} else{
			ss <- mxModel(model=model,
				mxMatrix('Zero', nksi, nksi, name='stateSpaceA'),
				mxMatrix('Zero', nksi, 1, name='stateSpaceX0'),
				mxMatrix('Iden', nksi, nksi, name='stateSpaceP0'),
				mxMatrix('Full', 1, 1, values=1, name='stateSpaceU'),
				mxExpectationStateSpace(A='stateSpaceA', B=model$expectation$KA, C=model$expectation$LX, D=model$expectation$TX, Q=model$expectation$PH, R=model$expectation$TD, x0='stateSpaceX0', P0='stateSpaceP0', u='stateSpaceU'))
			resDel <- mxKalmanScores(ss)
			res[,,1] <- resDel$xUpdated[-1,, drop=FALSE]
			res[,,2] <- apply(resDel$PUpdated[,,-1, drop=FALSE], 3, function(x){sqrt(diag(x))})

			# Kill the rows with too much missing data
			useCols <- dimnames(ss[[ss$expectation$C]])[[1]]
			rawData <- model$data$observed[ , useCols, drop=FALSE]
			res <- clearExcessivelyMissingRows(rawData, minManifests, res)
		}
	} else {
		stop('Unknown type argument to mxFactorScores')
	}
	dimnames(res) <- list(1:dim(res)[1], factorNames, c('Scores', 'StandardErrors'))
	return(res)
}

clearExcessivelyMissingRows <- function(rawData, minManifests, res) {
	numPresent <- apply(!is.na(rawData), 1, sum)
	if (any(numPresent < ncol(rawData)) && is.na(minManifests)) {
		requireMinManifests(which(numPresent < ncol(rawData))[1])
	}
	res[numPresent < minManifests, , 1] <- NA
	res[numPresent < minManifests, , 2] <- NA
	res
}

lisrelFactorScoreHelper <- function(model){
	lx <- mxEvalByName(model$expectation$LX, model, compute=TRUE)
	nksix <- dim(lx)
	nksi <- nksix[2]
	nx <- nksix[1]
	ksiMean <- mxEvalByName(model$expectation$KA, model, compute=TRUE)
	newKappa <- mxMatrix("Full", nksi, 1, values=ksiMean, free=TRUE, name="Score", labels=paste0("fscore", 1:nksi))
	scoreKappa <- mxAlgebraFromString(paste("Score -", model$expectation$KA), name="SKAPPA", dimnames=list(dimnames(lx)[[2]], 'one'))
	newExpect <- model$expectation
  newExpect$KA <- "SKAPPA"
  newExpect@.discreteCheckCount <- FALSE
	newWeight <- mxAlgebraFromString(paste0("log(det(", model$expectation$PH, ")) + ( (t(SKAPPA)) %&% ", model$expectation$PH, " ) + ", nksi, "*log(2*3.1415926535)"), name="weight")
	work <- mxModel(model=model, name=paste("FactorScores", model$name, sep=''), newKappa, scoreKappa, newExpect, newWeight)
	return(work)
}

createOppositeF <- function(Fmatrix){
	is.manifest <- as.logical(colSums(Fmatrix))
	mdim <- nrow(Fmatrix)
	tdim <- ncol(Fmatrix)
	ldim <- tdim - mdim
	tnam <- dimnames(Fmatrix)[[2]]
	lnam <- tnam[!is.manifest]
	OFmatrix <- matrix(0, nrow=ldim, ncol=tdim, dimnames=list(lnam, tnam))
	OFmatrix[lnam, lnam] <- diag(1, nrow=ldim)
	return(list(OF=OFmatrix, is.manifest=is.manifest))
}

ramFactorScoreHelper <- function(model){
	Fmat <- mxEvalByName(model$expectation$F, model, compute=TRUE)
	alldim <- dim(Fmat)
	tdim <- alldim[2]
	mdim <- alldim[1]
	ldim <- tdim - mdim
	OFmat <- createOppositeF(Fmat)
	fullMean <- mxEvalByName(model$expectation$M, model, compute=TRUE)
	scoreStart <- fullMean
	scoreStart[!OFmat$is.manifest] <- 0
	basVal <- fullMean
	basVal[OFmat$is.manifest] <- 0
	basNam <- paste0("Base", model$expectation$M)
	newMean <- mxMatrix("Full", 1, tdim, values=scoreStart, free=!OFmat$is.manifest, name="Score", labels=paste0("fscore", dimnames(Fmat)[[2]]))
	basMean <- mxMatrix("Full", 1, tdim, values=basVal, free=FALSE, name=basNam)
	scoreMean <- mxAlgebraFromString(paste("Score -", basNam), name="ScoreMinusM", dimnames=list('one', dimnames(Fmat)[[2]]))
	newExpect <- model$expectation
  newExpect@.discreteCheckCount <- FALSE
  newExpect$M <- "ScoreMinusM"
	oppF <- mxMatrix('Full', nrow=tdim-mdim, ncol=tdim, values=OFmat$OF, name='oppositeF')
	imat <- mxMatrix('Iden', tdim, tdim, name='IdentityMatrix')
	imaInv <- mxAlgebraFromString(paste("solve(IdentityMatrix - ", model$expectation$A, ")"), name='IdentityMinusAInverse')
	lcov <- mxAlgebraFromString(paste("oppositeF %*% IdentityMinusAInverse %*% ", model$expectation$S, " %*% t(IdentityMinusAInverse) %*% t(oppositeF)"), name='TheLatentRAMCovariance')
	newWeight <- mxAlgebraFromString(paste0("log(det(TheLatentRAMCovariance)) + ( (ScoreMinusM %*% t(oppositeF)) %&% TheLatentRAMCovariance ) + ", ldim, "*log(2*3.1415926535)"), name="weight")
	work <- mxModel(model=model, name=paste("FactorScores", model$name, sep=''), newMean, scoreMean, basMean, newExpect, oppF, imat, imaInv, lcov, newWeight)
	return(work)
}

RAMrfs <- function(model, res, minManifests) {
	i <- j <- 1
	manvars <- model@manifestVars
	latvars <- model@latentVars
	defvars <- findIntramodelDefVars(model)
	relevantDataCols <- c(manvars,defvars)
	dat <- model@data@observed
	I <- diag(length(manvars)+length(latvars))
	Ilat <- diag(length(latvars))
	while(i<=dim(res)[1]){
		continublockflag <- ifelse(i<dim(res)[1],TRUE,FALSE)
		manvars.curr <- manvars[ !is.na(dat[i,manvars]) ]
		while(continublockflag){
			#To include a subsequent row in the current block of rows,
			#we need to be sure that its missingness pattern is the same, and that if there are definition variables,
			#that their values in the subsequent row are equal to those in the previous rows:
			if(
				j<dim(res)[1] &&
				all(is.na(dat[j,relevantDataCols])==is.na(dat[(j+1),relevantDataCols])) &&
				( !length(defvars) || all(dat[j,defvars]==dat[(j+1),defvars]) )
			){j <- j+1}
			else{continublockflag <- FALSE}
		}
		unfilt <- solve(I-mxEvalByName("A",model,T,defvar.row=i))%*%mxEvalByName("S",model,T,defvar.row=i)%*%
			t(solve(I-mxEvalByName("A",model,T,defvar.row=i)))
		dimnames(unfilt) <- list(c(manvars,latvars),c(manvars,latvars)) #<--Necessary?
		latmeans <- matrix(1,ncol=1,nrow=(j-i+1)) %x% t(solve(Ilat-mxEvalByName("A",model,T,defvar.row=i)[latvars,latvars]) %*%
			matrix(mxEvalByName("M",model,T,defvar.row=i)[,latvars],ncol=1))
		missing <- is.na(dat[i,manvars])
		anyMissing <- any(missing)
		if (anyMissing && is.na(minManifests)) requireMinManifests(i)
		if (anyMissing && sum(!missing) < minManifests) {
				res[i:j,,1] <- NA
				res[i:j,,2] <- NA
		} else {
			if(all(missing)){
				res[i:j,,1] <- latmeans
				res[i:j,,2] <- matrix(1,ncol=1,nrow=(j-i+1)) %x% matrix(sqrt(diag(unfilt[latvars,latvars])),nrow=1)
			}
			else{
				obsmeans <- matrix(1,ncol=1,nrow=(j-i+1)) %x%
				matrix(mxGetExpected(model,"means",defvar.row=i)[,which(!is.na(dat[i,manvars]))],nrow=1)
				dat.curr <- as.matrix(dat[i:j,manvars.curr])
				if(i==j){dat.curr <- matrix(dat.curr,nrow=1)} #<--Annoying...
				res[i:j,,1] <- ( (dat.curr - obsmeans) %*%
						 (solve(unfilt[manvars.curr,manvars.curr])%*%unfilt[manvars.curr,latvars]) ) + latmeans
				indeterminateVariance <- unfilt[latvars,latvars] -
					(unfilt[latvars,manvars.curr]%*%solve(unfilt[manvars.curr,manvars.curr])%*%
					 unfilt[manvars.curr,latvars])
				res[i:j,,2] <- matrix(1,ncol=1,nrow=(j-i+1)) %x% matrix(sqrt(diag(indeterminateVariance)),nrow=1)
			}
		}
		i <- j+1
		j <- i
	}
	return(res)
}

findIntramodelDefVars <- function(model){
	if( length(model@runstate) && !length(model@runstate$defvars) ){return(NULL)}
	matlabs <- unlist(lapply(model@matrices,FUN=function(x){x@labels[!is.na(x@labels)]}))
	if( !("data." %in% substr(matlabs,1,5)) ){return(NULL)}
	else{
		defvars <- matlabs[which(substr(matlabs,1,5)=="data.")]
		defvars <- substr(defvars,6,nchar(defvars))
	}
	return( defvars )
}
