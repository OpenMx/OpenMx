#
#   Copyright 2007-2015 The OpenMx Project
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
# Author:
# Date: 2015-06-14
# Filename: MxFactorScores.R
# Purpose: Write a helper function for computing various type of factor scores
#------------------------------------------------------------------------------


mxFactorScores <- function(model, type=c('ML', 'WeightedML', 'Regression')){
	if(model$data$type!='raw'){
		stop("The 'model' arugment must have raw data.")
	}
	classExpect <- class(model$expectation)
	# TODO Add handling of multigroup models
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
		factorScoreHelperFUN <- lisrelFactorScoreHelper
	} else if(classExpect %in% "MxExpectationRAM"){
		fm <- mxEvalByName(model$expectation$F, model, compute=TRUE)
		nksix <- dim(fm)
		nksi <- nksix[2] - nksix[1]
		nx <- nksix[1]
		factorScoreHelperFUN <- ramFactorScoreHelper
	}
	nrows <- nrow(model$data$observed)
	res <- array(NA, c(nrows, nksi, 2))
	if(any(type %in% c('ML', 'WeightedML'))){
		model <- omxSetParameters(model, labels=names(omxGetParameters(model)), free=FALSE)
			work <- factorScoreHelperFUN(model)
		if(type[1]=='WeightedML'){
			wup <- mxModel(model="Container", work,
				mxAlgebraFromString(paste(work@name, ".weight + ", work@name, ".fitfunction", sep=""), name="wtf"),
				mxFitFunctionAlgebra("wtf")
			)
			work <- wup
		}
		work@data <- NULL
		for(i in 1:nrows){
			if(type[1]=='ML'){
				fit <- mxModel(model=work, name=paste0(work@name, i, "Of", nrows), mxData(model$data$observed[i,,drop=FALSE], 'raw'))
			} else if(type[1]=='WeightedML'){
				work@submodels[[1]]@data <- mxData(model$data$observed[i,,drop=FALSE], 'raw')
				fit <- mxModel(model=work, name=paste0(work@name, i, "Of", nrows))
			}
			fit <- mxRun(fit, silent=as.logical((i-1)%%100), suppressWarnings=TRUE)
			res[i,,1] <- omxGetParameters(fit) #params
			res[i,,2] <- fit$output$standardErrors #SEs
		}
	} else if(type=='Regression'){
		if(!single.na(model$expectation$thresholds)){
			stop('Regression factor scores cannot be computed when there are thresholds (ordinal data).')
		}
		if(!(classExpect %in% "MxExpectationLISREL")){
			stop('Regression factor scores are only possible for LISREL expectations.')
		}
		ss <- mxModel(model=model,
			mxMatrix('Zero', nksi, nksi, name='stateSpaceA'),
			mxMatrix('Zero', nksi, nx, name='stateSpaceB'),
			mxMatrix('Iden', nx, nx, name='stateSpaceD'),
			mxMatrix('Iden', nksi, nksi, name='stateSpaceP0'),
			mxExpectationStateSpace(A='stateSpaceA', B='stateSpaceB', C=model$expectation$LX, D='stateSpaceD', Q=model$expectation$PH, R=model$expectation$TD, x0=model$expectation$KA, P0='stateSpaceP0', u=model$expectation$TX))
		resDel <- mxKalmanScores(ss)
		res[,,1] <- resDel$xUpdated[-1,, drop=FALSE]
		res[,,2] <- apply(resDel$PUpdated[,,-1, drop=FALSE], 3, function(x){sqrt(diag(x))})
	} else {
		stop('Unknown type argument to mxFactorScores')
	}
	return(res)
}

lisrelFactorScoreHelper <- function(model){
	lx <- mxEvalByName(model$expectation$LX, model, compute=TRUE)
	nksix <- dim(lx)
	nksi <- nksix[2]
	nx <- nksix[1]
	ksiMean <- mxEvalByName(model$expectation$KA, model, compute=TRUE)
	newKappa <- mxMatrix("Full", nksi, 1, values=ksiMean, free=TRUE, name="Score", labels=paste0("fscore", 1:nksi))
	scoreKappa <- mxAlgebraFromString(paste("Score -", model$expectation$KA), name="SKAPPA", dimnames=list(dimnames(lx)[[2]], 'one'))
	newExpect <- mxExpectationLISREL(LX=model$expectation$LX, PH=model$expectation$PH, TD=model$expectation$TD, TX=model$expectation$TX, KA="SKAPPA", thresholds=model$expectation$thresholds)
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
	newMean <- mxMatrix("Full", 1, tdim, values=scoreStart, free=!OFmat$is.manifest, name="Score", labels=paste0("fscore", dimnames(Fmat)[[2]]))
	scoreMean <- mxAlgebraFromString(paste("Score -", model$expectation$M), name="ScoreMinusM", dimnames=list('one', dimnames(Fmat)[[2]]))
	newExpect <- mxExpectationRAM(A=model$expectation$A, S=model$expectation$S, F=model$expectation$F, M="ScoreMinusM", thresholds=model$expectation$thresholds)
	oppF <- mxMatrix('Full', nrow=tdim-mdim, ncol=tdim, values=OFmat$OF, name='oppositeF')
	imat <- mxMatrix('Iden', tdim, tdim, name='IdentityMatrix')
	imaInv <- mxAlgebraFromString(paste("solve(IdentityMatrix - ", model$expectation$A, ")"), name='IdentityMinusAInverse')
	lcov <- mxAlgebraFromString(paste("oppositeF %*% IdentityMinusAInverse %*% ", model$expectation$S, " %*% t(IdentityMinusAInverse) %*% t(oppositeF)"), name='TheLatentRAMCovariance')
	newWeight <- mxAlgebraFromString(paste0("log(det(TheLatentRAMCovariance)) + ( (ScoreMinusM %*% t(oppositeF)) %&% TheLatentRAMCovariance ) + ", ldim, "*log(2*3.1415926535)"), name="weight")
	work <- mxModel(model=model, name=paste("FactorScores", model$name, sep=''), newMean, scoreMean, newExpect, oppF, imat, imaInv, lcov, newWeight)
	return(work)
}

