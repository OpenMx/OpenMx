#   Copyright 2007-2018 The OpenMx Project
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

#mxMMI <- function(reference=character(0), models=list(), covariances=NULL, include=c("onlyFree","all"), AICc=FALSE, confidLvl=0.95){
omxAkaikeWeights <- function(models=list(), type=match.arg("AIC","AICc"), confidLvl=0.95){	
	#Input checking:
	type <- match.arg(type,c("AIC","AICc"))
	if(length(models)<2 || !all(sapply(models,is,"MxModel"))){
		stop("argument 'models' must be a list of at least 2 MxModel objects")
	}
	if(confidLvl>=1 || confidLvl<=0){stop("argument 'confidLvl' must be a proportion strictly between 0 and 1")}

	modelNames <- rep(NA_character_, length(models))
	AICs <- rep(NA_real_, length(models))
	isGREML <- NULL
	gfxf <- NULL
	#Loop thru models to do safety & validity checks and calculate AIC for each:
	for(i in 1:length(models)){
		warnModelCreatedByOldVersion(models[[i]])
		if(.hasSlot(models[[i]],".wasRun") && !models[[i]]@.wasRun){
			stop(paste("MxModel",omxQuotes(models[[i]]@name),"has not yet been run"))
		}
		if(!length(models[[i]]@output)){
			stop(paste("MxModel",omxQuotes(models[[i]]@name),"has an empty 'output' slot; has it been run yet?"))
		}
		if(.hasSlot(models[[i]],".modifiedSinceRun") && models[[i]]@.modifiedSinceRun){
			warning(paste("MxModel",omxQuotes(models[[i]]@name),"has been modified since it was run"))
		}
		if(i>1 && (models[[i]]@name %in% modelNames[1:(i-1)])){
			stop("each MxModel object in 'models' must be uniquely identified by the character string in its 'name' slot")
		}
		if( length(models[[i]]@output$fitUnits) && models[[i]]@output$fitUnits!="-2lnL" ){
			stop(paste("MxModel",omxQuotes(models[[i]]@name),"does not have -2lnL fit units"))
		}
		currsumm <- summary(models[[i]])
		if(is(models[[i]]@fitfunction,"MxFitFunctionGREML") ){
			isGREML <- c(isGREML,TRUE)
			gfxf <- c(gfxf, paste(currsumm$GREMLfixeff$name,collapse=","))
			if(length(unique(isGREML))>1){stop("some but not all of 'models' use MxFitFunctionGREML")}
		} 
		else{
			isGREML <- c(isGREML,FALSE)
		}
		modelNames[i] <- models[[i]]@name
		if(type=="AIC"){
			AICs[i] <- currsumm$informationCriteria["AIC:","par"]
		}	
		if(type=="AICc"){
			AICs[i] <- currsumm$informationCriteria["AIC:","sample"]
		}
		if(!is.finite(AICs[i])){
			stop(paste("MxModel",omxQuotes(models[[i]]@name),"has a non-finite AIC"))
		}
	}
	
	if(length(gfxf) && length(unique(gfxf))>1){
		#As with mxCompare(), this is a warning, not an error:
		warning("not all of 'models' have matching names for their fixed effects; results may be invalid")
	}
	
	deltas <- AICs - min(AICs)
	#Akaike weights:
	w <- exp(-deltas/2) / sum(exp(-deltas/2))
	
	#Output table
	tabl1 <- data.frame(model=modelNames,AIC=AICs,delta=deltas,AkaikeWeight=w,inConfidenceSet=rep("",length(models)),stringsAsFactors=F)
	if(type=="AICc"){colnames(tabl1)[2] <- "AICc"}
	tabl1 <- tabl1[order(tabl1[,2]),]
	cump <- 0
	for(i in 1:length(models)){
		tabl1$inConfidenceSet[i] <- "*"
		cump <- cump + tabl1$AkaikeWeight[i]
		if(cump >= confidLvl){break}
	}
	attr(tabl1,"unsortedModelNames") <- modelNames
	return(tabl1)
}


mxModelavgEstimates <- 	function(
	reference=character(0), models=list(), include=c("onlyFree","all"), asBlock=FALSE, covariances=list(), type=c("AIC","AICc"), 
	confidLvl=0.95)
{
	include <- match.arg(include,c("onlyFree","all"))
	type <- match.arg(type,c("AIC","AICc"))
	tabl1 <- omxAkaikeWeights(models=models, type=type, confidLvl=confidLvl)
	names(models) <- attr(tabl1,"unsortedModelNames")
	models <- models[tabl1$"model"]
	if(length(covariances)){stop("Not Yet Implemented")}
	if(asBlock){
		stop("Not Yet Implemented")
	} 
	
	refdims <- matrix(NA_real_,nrow=length(reference),ncol=2,dimnames=list(reference,c("nrow","ncol")))
	modreflist <- vector("list",length(models))
	names(modreflist) <- names(models)
	for(i in 1:length(modreflist)){
		currmod <- models[[i]]
		modreflist[[i]] <- vector("list",length(reference))
		names(modreflist[[i]]) <- reference
		for(j in 1:length(reference)){
			xx <- try(mxEvalByName(name=reference[j],model=currmod,compute=T))
			if(is(xx,"try-error") || !length(xx)){next} #<--Maybe throw warning
			if(all(is.na(refdims[j,]))){refdims[j,] <- c(nrow(xx),ncol(xx))}
			else if( !(refdims[j,1]==nrow(xx) && refdims[j,2]==ncol(xx)) ){
				stop(paste("reference ",omxQuotes(reference[j])," is ",nrow(xx),"x",ncol(xx)," in MxModel ",omxQuotes(currmod@name)," but is ",refdims[j,1],"x",refdims[j,2]," in at least one other MxModel",sep=""))
			}
			modreflist[[i]][[j]] <- xx
		}
	}
	reflengths <- refdims[,1]*refdims[,2]
	
	longlabels <- NULL
	for(i in 1:nrow(refdims)){
		if(refdims[i,1]==1 && refdims[i,2]==1){
			longlabels <- c(longlabels, reference[i])
		}
		else{
			for(j in 1:refdims[i,2]){
				for(k in 1:refdims[i,1]){
					longlabels <- c(longlabels, paste(reference[i],"[",k,",",j,"]",sep=""))
				}
			}
		}
	}
	
	
	thetamtx <- matrix(NA_real_,nrow=length(longlabels),ncol=nrow(tabl1),dimnames=list(longlabels,names(models)))
	wivmtx <- matrix(NA_real_,nrow=length(longlabels),ncol=nrow(tabl1),dimnames=list(longlabels,names(models)))
	# if(asBlock){
	# 	stop("Not Yet Implemented")
	# } 
	#else{
	for(i in 1:length(models)){
		currmod <- models[[i]]
		currmodparams <- list(fre=omxGetParameters(currmod,free=TRUE),fix=omxGetParameters(currmod,free=FALSE))
		#TODO: check here for user-supplied cov matrix
		if(!is.na(currmod@output$infoDefinite) && currmod@output$infoDefinite){
			currmodcovm <- 2*chol2inv(chol(currmod@output$hessian))
			dimnames(currmodcovm) <- dimnames(currmod@output$hessian)
		}
		else{currmodcovm <- 2*solve(currmod@output$hessian)}
		rownumcurr <- 1
		for(j in 1:length(reference)){
			if(reference[j] %in% names(currmodparams$fre)){
				thetamtx[rownumcurr,i] <- currmodparams$fre[reference[j]]
				wivmtx[rownumcurr,i] <- currmodcovm[reference[j],reference[j]]
				rownumcurr <- rownumcurr + 1
				next
			} 
			else if(reference[j] %in% names(currmodparams$fix)){
				if(include=="onlyFree"){
					rownumcurr <- rownumcurr + 1
					next
				}
				thetamtx[rownumcurr,i] <- currmodparams$fix[reference[j]]
				wivmtx[rownumcurr,i] <- 0
				rownumcurr <- rownumcurr + 1
				next
			}
			else{ #i.e., if it's not a parameter label
				x <- modreflist[[i]][[j]]
				if(!length(x)){
					rownumcurr <- rownumcurr + reflengths[j]
					next
				}
				#thetamtx[rownumcurr:(rownumcurr+reflengths[j]-1),i] <- x
				xv <- mxSE(x=reference[j],model=currmod,details=T,cov=currmodcovm,forceName=T)$Cov
				for(k in 1:nrow(xv)){
					if(xv[k] < .Machine$double.eps && include=="onlyFree"){next}
					thetamtx[rownumcurr+(k-1),i] <- x[k]
					wivmtx[rownumcurr+(k-1),i] <- xv[k,k]
				}
				rownumcurr <- rownumcurr + reflengths[j]
			}
		}
	}
	#return(list(thetamtx,wivmtx))
	
	tabl2 <- matrix(NA_real_,ncol=2,nrow=nrow(thetamtx),dimnames=list(longlabels,c("Estimate","SE")))
	for(i in 1:nrow(tabl2)){
		thetacurr <- thetamtx[i,]
		wivcurr <- wivmtx[i,]
		wcurr <- tabl1$AkaikeWeight
		if(any(is.na(thetacurr))){
			thetacurr <- thetacurr[which(!is.na(thetamtx[i,]))]
			wivcurr <- wivcurr[which(!is.na(thetamtx[i,]))]
			wcurr <- wcurr[which(!is.na(thetamtx[i,]))]
			wcurr <- wcurr/sum(wcurr)
		}
		tabl2[i,1] <- t(wcurr) %*% thetacurr
		tabl2[i,2] <- sqrt( t(wcurr) %*% (wivcurr + (tabl2[i,1]-thetacurr)^2) )
	}
	return(list(tabl2,thetamtx,wivmtx,tabl1))
}



