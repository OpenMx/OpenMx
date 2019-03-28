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

#mxMMI <- function(reference=character(0), models=list(), covariances=NULL, include=c("onlyFree","all"), AICc=FALSE, conf.level=0.95){
omxAkaikeWeights <- function(models=list(), type=c("AIC","AICc"), conf.level=0.95){	
	#Input checking:
	type <- match.arg(type,c("AIC","AICc"))
	if(length(models)<2 || !all(sapply(models,is,"MxModel"))){
		stop("argument 'models' must be a list of at least 2 MxModel objects")
	}
	conf.level <- conf.level[1]
	if(conf.level>=1 || conf.level<=0){stop("argument 'conf.level' must be a proportion strictly between 0 and 1")}
	
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
		if(cump >= conf.level){break}
	}
	attr(tabl1,"unsortedModelNames") <- modelNames
	return(tabl1)
}


mxModelAverage <- function(
	reference=character(0), models=list(), include=c("onlyFree","all"), SE=NULL, refAsBlock=FALSE, covariances=list(), type=c("AIC","AICc"), 
	conf.level=0.95)
{
	#Input checking:
	include <- match.arg(include,c("onlyFree","all"))
	type <- match.arg(type,c("AIC","AICc"))
	if(!length(SE)){
		SE <- ifelse(include=="onlyFree",TRUE,FALSE)
	}
	#The list of models is thoroughly checked as part of this function call:
	tabl1 <- omxAkaikeWeights(models=models, type=type, conf.level=conf.level)
	#If 'reference' is empty, return tabl1, with a warning:
	if(!length(reference)){
		warning("argument 'reference' is NULL; vector of character strings was expected")
		return(tabl1)
	}
	names(models) <- attr(tabl1,"unsortedModelNames")
	models <- models[tabl1$"model"] #<--Rearrange 'models' in order from best to worst AIC.
	#if covariances is non-empty, it must either be named, or be of the same length as 'models':
	if(length(covariances)){
		if(!length(names(covariances))){
			if(length(covariances) != length(models)){
				stop("non-empty value for argument 'covariances' must be a named list, or a list of the same length as argument 'models'")
			}
			names(covariances) <- attr(tabl1,"unsortedModelNames")
			covariances <- covariances[tabl1$model]
		}
		else{
			if(length(unique(names(covariances))) < length(names(covariances))){
				stop("if a non-empty value for argument 'covariances' has element names, then those names must uniquely identify each element")
			}
			if( !all(names(covariances) %in% names(models)) ){
				#TODO more helpful warning message:
				warning(paste("an element of argument 'covariances' has a name that does not match the name of any element of 'models'"))
			}
		}
	}
	
	#Function to be used to identify fixed algebra elements:
	sefun <- function(x = NULL, model, alg){
		if(is.null(x)){x <- omxGetParameters(model)}
		param <- omxGetParameters(model)
		paramNames <- names(param)
		model <- omxSetParameters(model, values=x, labels=paramNames, free=TRUE)
		out <- try(mxEvalByName(alg, model, compute=TRUE),silent=T)
		if(is(out,"try-error")){
			model <- mxModel(model,mxAlgebraFromString(algString=alg,name="onTheFlyAlgebra2"))
			out <- mxEvalByName(name="onTheFlyAlgebra2",model=model,compute=T)
		}
		return(out)
	}
	
	#Now, we need to find out how many elements long each reference is, and make sure that each reference refers to a matrix of the same
	#dimensions in all models in which it can be evaluated.  We'll store the results of mxEvalByName() and later put them into a matrix:
	refdims <- matrix(NA_real_,nrow=length(reference),ncol=2,dimnames=list(reference,c("nrow","ncol")))
	modreflist <- vector("list",length(models))
	names(modreflist) <- names(models)
	for(i in 1:length(modreflist)){
		currmod <- models[[i]]
		modreflist[[i]] <- vector("list",length(reference))
		names(modreflist[[i]]) <- reference
		for(j in 1:length(reference)){
			xx <- try(mxEvalByName(name=reference[j],model=currmod,compute=T),silent=TRUE)
			if(is(xx,"try-error")){
				currmod <- mxModel(currmod,mxAlgebraFromString(algString=reference[j],name="onTheFlyAlgebra2"))
				xx <- try(mxEvalByName(name="onTheFlyAlgebra2",model=currmod,compute=T),silent=TRUE)
			}
			if(is(xx,"try-error") || !length(xx)){
				if(refAsBlock){stop(paste("reference",omxQuotes(reference[j]),"could not be evaluated in model",omxQuotes(currmod@name)))}
				else{
					warning(paste("reference",omxQuotes(reference[j]),"could not be evaluated in model",omxQuotes(currmod@name)))
					next
				}
			}
			if(all(is.na(refdims[j,]))){refdims[j,] <- c(nrow(xx),ncol(xx))}
			else if( !(refdims[j,1]==nrow(xx) && refdims[j,2]==ncol(xx)) ){
				stop(paste("reference ",omxQuotes(reference[j])," is ",nrow(xx),"x",ncol(xx)," in MxModel ",omxQuotes(currmod@name)," but is ",refdims[j,1],"x",refdims[j,2]," in at least one other MxModel",sep=""))
			}
			modreflist[[i]][[j]] <- xx
		}
	}
	reflengths <- refdims[,1]*refdims[,2]
	
	#Make labels for elements of matrices and algebras:
	longlabels <- NULL
	for(i in 1:nrow(refdims)){
		#Possible TODO--instead of this error, drop the bad reference and carry on:
		if(any(is.na(refdims[i,]))){
			stop(paste("reference ",omxQuotes(rownames(refdims)[i])," could not be evaluated in any model",sep=""))
		}
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
	if(refAsBlock){
		refcovlist <- vector("list",length(models))
		for(i in 1:length(models)){
			currmod <- models[[i]]
			currmodparams <- list(fre=omxGetParameters(currmod,free=TRUE),fix=omxGetParameters(currmod,free=FALSE))
			covmAvailable <- FALSE
			if(SE){
				if(currmod@name %in% names(covariances)){
					currmodcovm <- covariances[[currmod@name]]
					#Sanity checks on user-supplied covariance matrices:
					if(!isSymmetric(currmodcovm)){ #<--Returns FALSE for non-square matrices
						stop(paste("covariance matrix for model",omxQuotes(currmod@name),"is not symmetric"))
					}
					if( !length(rownames(currmodcovm)) || !length(colnames(currmodcovm)) ){
						stop(paste("covariance matrix for model",omxQuotes(currmod@name),"does not have both row names and column names"))
					}
					if( !all(rownames(currmodcovm) == colnames(currmodcovm)) ){
						stop(paste("the row names and column names of covariance matrix for model",omxQuotes(currmod@name),"are not equal"))
					}
					if( !all(rownames(currmodcovm) %in% names(currmodparams$fre)) ||  
							!all(names(currmodparams$fre) %in% rownames(currmodcovm)) ){
						stop(paste("the dimnames of the covariance matrix for model",omxQuotes(currmod@name),"do not match the labels of the model's free parameters"))
					}
					covmAvailable <- TRUE
				}
				else{
					currmodcovm <- try(vcov(currmod))
					if(!is(currmodcovm,"try-error")){covmAvailable <- TRUE}
					# if(!is.na(currmod@output$infoDefinite) && currmod@output$infoDefinite){
					# 	currmodcovm <- 2*chol2inv(chol(currmod@output$hessian))
					# 	dimnames(currmodcovm) <- dimnames(currmod@output$hessian)
					# 	covmAvailable <- TRUE
					# }
					# else{
					# 	currmodcovm <- try(solve(currmod@output$hessian/2))
					# 	if(!is(currmodcovm,"try-error")){covmAvailable <- TRUE}
					# }
					# if(imxHasConstraint(currmod) && covmAvailable){
					# 	warning(paste("due to presence of MxConstraints in model",omxQuotes(currmod@name),"sampling covariance matrix and standard errors for model-average point estimates may not be valid"))
					# }
				}
			}
			currmod <- mxModel(
				currmod,
				mxAlgebraFromString(algString=paste("rbind(",paste(longlabels,collapse=","),")",sep=""),name="onTheFlyAlgebra2"))
			thetamtx[,i] <- mxEvalByName(name="onTheFlyAlgebra2",model=currmod,compute=T)[,1]
			if(!covmAvailable && SE){
				#Not sure what else to do in this case:
				warning(paste("model",omxQuotes(currmod@name),"has no covariance matrix for its free parameters; continuing with argument 'SE' coerced to FALSE"))
				SE <- FALSE
			}
			if(covmAvailable && SE){
				refcov <- mxSE(x="onTheFlyAlgebra2",model=currmod,details=T,cov=currmodcovm,silent=T)$Cov
				#TODO more informative error message:
				if(include=="onlyFree" && any(diag(refcov) <= .Machine$double.eps)){
					stop("when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
				}
				refcovlist[[i]] <- refcov
			}
			else if(include=="onlyFree"){
				jac <- numDeriv::jacobian(func=sefun,x=omxGetParameters(currmod),model=currmod,alg="onTheFlyAlgebra2")
				if( any(apply(X=jac,MARGIN=1,FUN=function(x){all(x==0)})) ){
					#TODO more informative error message:
					stop("when refAsBlock=TRUE and include='onlyFree', no references may be fixed in any model")
				}
			}
		}
		tabl2 <- matrix(0.0,ncol=2,nrow=nrow(thetamtx),dimnames=list(longlabels,c("Estimate","SE")))
		tabl2[,1] <- thetamtx%*%matrix(tabl1$AkaikeWeight,ncol=1)
		if(SE){
			uncondcovm <- matrix(0.0,nrow=nrow(thetamtx),ncol=nrow(thetamtx),dimnames=list(longlabels,longlabels))
			w <- tabl1$AkaikeWeight
			for(i in 1:length(refcovlist)){
				bw <- as.vector(thetamtx[,i]-tabl2[,1])
				uncondcovm <- uncondcovm + tabl1$AkaikeWeight[i]*(refcovlist[[i]] + outer(bw,bw))
			}
			tabl2[,2] <- sqrt(diag(uncondcovm))
			outlist <- list(tabl2,thetamtx,uncondcovm,tabl1)
			names(outlist) <- c("Model-Average Estimates","Model-wise Estimates","Joint Covariance Matrix","Akaike-Weights Table")
			return(outlist)
		}
		else{
			tabl2 <- matrix(tabl2[,1],nrow=nrow(tabl2),ncol=1,dimnames=list(longlabels,"Estimate"))
			outlist <- list(tabl2,thetamtx,NULL,tabl1)
			names(outlist) <- c("Model-Average Estimates","Model-wise Estimates","Joint Covariance Matrix","Akaike-Weights Table")
			return(outlist)
		}
	}
	else{
		wivmtx <- matrix(NA_real_,nrow=length(longlabels),ncol=nrow(tabl1),dimnames=list(longlabels,names(models)))
		for(i in 1:length(models)){
			currmod <- models[[i]]
			currmodparams <- list(fre=omxGetParameters(currmod,free=TRUE),fix=omxGetParameters(currmod,free=FALSE))
			covmAvailable <- FALSE
			if(currmod@name %in% names(covariances)){
				currmodcovm <- covariances[[currmod@name]]
				#Sanity checks on user-supplied covariance matrices:
				if(!isSymmetric(currmodcovm)){ #<--Returns FALSE for non-square matrices
					stop(paste("covariance matrix for model",omxQuotes(currmod@name),"is not symmetric"))
				}
				if( !length(rownames(currmodcovm)) || !length(colnames(currmodcovm)) ){
					stop(paste("covariance matrix for model",omxQuotes(currmod@name),"does not have both row names and column names"))
				}
				if( !all(rownames(currmodcovm) == colnames(currmodcovm)) ){
					stop(paste("the row names and column names of covariance matrix for model",omxQuotes(currmod@name),"are not equal"))
				}
				if( !all(rownames(currmodcovm) %in% names(currmodparams$fre)) ||  
						!all(names(currmodparams$fre) %in% rownames(currmodcovm)) ){
					stop(paste("the dimnames of the covariance matrix for model",omxQuotes(currmod@name),"do not match the labels of the model's free parameters"))
				}
				covmAvailable <- TRUE
			}
			else{
				currmodcovm <- try(vcov(currmod))
				if(!is(currmodcovm,"try-error")){covmAvailable <- TRUE}
				# if(!is.na(currmod@output$infoDefinite) && currmod@output$infoDefinite){
				# 	currmodcovm <- 2*chol2inv(chol(currmod@output$hessian))
				# 	dimnames(currmodcovm) <- dimnames(currmod@output$hessian)
				# 	covmAvailable <- TRUE
				# }
				# else{
				# 	currmodcovm <- try(solve(currmod@output$hessian/2))
				# 	if(!is(currmodcovm,"try-error")){covmAvailable <- TRUE}
				# }
				# if(imxHasConstraint(currmod) && covmAvailable){
				# 	if(SE){
				# 		warning(paste("due to presence of MxConstraints in model",omxQuotes(currmod@name),"its model-wise sampling variances, as well as the standard errors for model-average point estimates, may not be valid"))
				# 	}
				# 	else{
				# 		warning(paste("due to presence of MxConstraints in model",omxQuotes(currmod@name),"its model-wise sampling variances may not be valid"))
				# 	}
				# }
			}
			if(!covmAvailable && SE){
				#We don't want the subset of models contributing to the model-average point estimates to be different from the subset of models 
				#contributing to their standard errors:
				warning(paste("model",omxQuotes(currmod@name),"has no covariance matrix for its free parameters; continuing with argument 'SE' coerced to FALSE"))
				SE <- FALSE
			}
			rownumcurr <- 1
			for(j in 1:length(reference)){
				if(reference[j] %in% names(currmodparams$fre)){
					thetamtx[rownumcurr,i] <- currmodparams$fre[reference[j]]
					if(covmAvailable){wivmtx[rownumcurr,i] <- currmodcovm[reference[j],reference[j]]}
					rownumcurr <- rownumcurr + 1
					next
				} 
				else if(reference[j] %in% names(currmodparams$fix)){
					if(include=="onlyFree"){
						rownumcurr <- rownumcurr + 1
						next
					}
					thetamtx[rownumcurr,i] <- currmodparams$fix[reference[j]]
					if(covmAvailable){wivmtx[rownumcurr,i] <- 0}
					rownumcurr <- rownumcurr + 1
					next
				}
				else{ #i.e., if it's not a parameter label
					x <- modreflist[[i]][[j]]
					if(!length(x)){
						rownumcurr <- rownumcurr + reflengths[j]
						next
					}
					if(covmAvailable){
						xv <- try(mxSE(x=reference[j],model=currmod,details=T,cov=currmodcovm,forceName=T,silent=T),silent=T)
						if(is(xv,"try-error")){
							currmod <- mxModel(currmod,mxAlgebraFromString(algString=reference[j],name="onTheFlyAlgebra2"))
							xv <- mxSE(x="onTheFlyAlgebra2",model=currmod,details=T,cov=currmodcovm,silent=T)
						}
						xv <- xv$Cov
						for(k in 1:nrow(xv)){
							if(include=="onlyFree" && xv[k,k] < .Machine$double.eps){next}
							thetamtx[rownumcurr+(k-1),i] <- x[k]
							wivmtx[rownumcurr+(k-1),i] <- xv[k,k]
						}
					}
					else{
						#We only care about detecting fixed algebra elements if include=="onlyFree":
						if(include=="all"){jac <- matrix(1,ncol=1,nrow=reflengths[j])}
						else{jac <- numDeriv::jacobian(func=sefun,x=omxGetParameters(currmod),model=currmod,alg=reference[j])}
						for(k in 1:reflengths[j]){
							#A row of the Jacobian matrix can only be all zeroes if include="onlyFree":
							if(all(jac[k,]==0)){next}
							thetamtx[rownumcurr+(k-1),i] <- x[k]
						}
					}
					rownumcurr <- rownumcurr + reflengths[j]
				}
			}
		}
		if(SE){
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
		}
		else{
			tabl2 <- matrix(NA_real_,ncol=1,nrow=nrow(thetamtx),dimnames=list(longlabels,c("Estimate")))
			for(i in 1:nrow(tabl2)){
				thetacurr <- thetamtx[i,]
				wcurr <- tabl1$AkaikeWeight
				if(any(is.na(thetacurr))){
					thetacurr <- thetacurr[which(!is.na(thetamtx[i,]))]
					wcurr <- wcurr[which(!is.na(thetamtx[i,]))]
					wcurr <- wcurr/sum(wcurr)
				}
				tabl2[i,1] <- t(wcurr) %*% thetacurr
			}
		}
		outlist <- list(tabl2,thetamtx,wivmtx,tabl1)
		names(outlist) <- c("Model-Average Estimates","Model-wise Estimates","Model-wise Sampling Variances","Akaike-Weights Table")
		return(outlist)
	}
}



