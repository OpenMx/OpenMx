omxTransformModelPPML <- function(model) {
	###############
	#CHECK SECTION#
	###############
	#checks are ordered from fast to slow
	#is the model a RAM model?
	objective <- model$objective
	if(is.null(objective) || !is(objective, "MxRAMObjective")) {
		return(model)
	}
	Aname <- objective$A
	Sname <- objective$S
	Fname <- objective$F
	Amatrix <- model[[Aname]]
	Smatrix <- model[[Sname]]
	Fmatrix <- model[[Fname]]
	if(!is(Amatrix, "MxMatrix") || !is(Smatrix, "MxMatrix") || !is(Fmatrix, "MxMatrix")) {
		return(model)
	}
	#transformation seems not really valuable in the situation of cov data
	#and is not written in this way to this point
	if(model@data@type != "raw") {
		return(model)
	}
	#are all regression loadings fixed?
	if(any(Amatrix@free)) {
		return(model)	
	}
	#do all errors have the same label?
	clabels <- array(0,0)
	fakeLatents <- array(0,0)
	manHasVar <- matrix(FALSE,length(model@manifestVars),1)
	rownames(manHasVar) <- model@manifestVars
	#gather all labels from the direct errors
	for(manifestVar in model@manifestVars) {
		if (Smatrix@values[manifestVar,manifestVar] !=0) {
			manHasVar[manifestVar,1] <- TRUE
			clabels <- append(clabels,model$S@labels[manifestVar,manifestVar])
		}		
	}
	#get all fakeLatens and their label
	for(latentVar in model@latentVars){
		#does the latentVar only have one pointer? => fakeLatent
		tmp <- which(Amatrix@values[model@manifestVars,latentVar] != 0)
		if (length(tmp) == 1) {
			#the regression weight has to be one
			if ((Amatrix@values[tmp,latentVar] == 1) && !manHasVar[tmp,1]) {
				fakeLatents <- append(fakeLatents,latentVar)
			} else {
				return(model)
			}
		}
	}
	for(fakeLatent in fakeLatents){
		clabels <- append(clabels, Smatrix@labels[fakeLatent,fakeLatent])
	}
	#not every manifest error variance has the same label
	if (length(unique(clabels)) != 1) {
		return(model)
	}
	#no latent is allowed to point at another latent
	if (any(Amatrix@values[model@latentVars,model@latentVars] != 0)) {
		return(model)
	}

	###################
	#TRANSFORM SECTION#
	###################
	#get loadings from latents to manifest from A,S,F
	#transorm A to E
	E <- solve(diag(nrow(Amatrix)) - Amatrix@values)
	#select only these columns of A which are real latents
	realLatents <- model@latentVars[is.na(pmatch(model@latentVars, fakeLatents))]
	#get loadings matrix:= The part of A which goes from the latents to the manifests
	lambda <- as.matrix(E[model@manifestVars, realLatents])
	qrDecom <- qr(lambda)
	#check if loadings matrix is of full rank
	k <- length(realLatents)
	if (qrDecom$rank != dim(lambda)[2] || k != qrDecom$rank){
		return(model)
	}
	#last check passed

	#ensure that every free parameter has a unique label
	pair <- omxNameAnonymousParameters(model)
	model <- pair[[1]]
	oldLabels <- pair[[2]]
	
	

	#calcualte upper triangle matrix
	#orthotogonal
	Q <- t(qr.Q(qrDecom, complete = TRUE))
	#transform loadings matrixn
	lambda <- Q %*% lambda
	#set all rows from k+1 to zero
	lambda[(k + 1) : nrow(lambda), ] <- 0
	#insert new loadings matrix in A
	Amatrix@values[model@manifestVars,realLatents] <- lambda
	model[[Aname]] <- Amatrix
	#tranform data or cov
	if(model$data@type == "raw"){
		model$data@observed <- as.matrix(model@data@observed) %*% t(Q)
		colnames(model$data@observed) <- model@manifestVars
	}
	
	###############
	#SPLIT SECTION#
	###############
	#leftmodel
	#calculate A realLatents indices for left model
	#first k manifest vars
	selectManifests <- model@manifestVars[1:k]
	#all "REAL" latents + fake latents which point into the selectManifests
	#remember that we checked that each fakeLatent column in A only has 1 non zero entry
	#which is one
	indices <- colSums(Amatrix@values[selectManifests, fakeLatents, drop = FALSE])

	selectFake <- fakeLatents[which(indices == 1)]
	selectLatents <- append(realLatents, selectFake)	
	leftmodel <- mxRename(model, 'leftmodel')	
	leftmodel <- selectSubModelFData(leftmodel, selectLatents, selectManifests)

	
	#rightmodel
	selectLatents <- model@latentVars[is.na(pmatch(model@latentVars, selectLatents))]
	selectManifests <- model@manifestVars[is.na(pmatch(model@manifestVars, selectManifests))]
	rightmodel <- mxRename(model, 'rightmodel')	
	rightmodel <- selectSubModelFData(rightmodel, selectLatents, selectManifests)

	
	#reunite the two submodel
	result <-  mxModel('PPMLModel', leftmodel, rightmodel)
	if (model$data@type == 'raw') {
		result$data <- model$data
	}
	modelnames <- c("leftmodel", "rightmodel")
	objectives <- paste(modelnames, "objective", sep = ".")
	objectives <- paste(objectives, collapse = " + ")
	expression <- paste("mxAlgebra(", objectives, ", name = 'TotObj')", sep = "")
	algebra <- eval(parse(text=expression))
	objective <- mxAlgebraObjective("TotObj")
	result <- mxModel(result, objective, algebra)
	return(result)
}

jlRun <- function(model){
	model <- omxTransformModelPPML(model)
	return(mxRun(model))
}
	

#@author:	Julian Karch
#		jk3nq@virginia.edu
#@idea>		Timo von Oertzen
#		timo@virginia.edu
selectSubModelFData <- function(model, selectLatents, selectManifests) {
	Aindices <- append(selectManifests, selectLatents)
	#build the two new models
	submodel <- model
	submodel@manifestVars <- selectManifests
	submodel@latentVars <- selectLatents
	submodel$A <- model$A[Aindices,Aindices]
	dimnames(submodel$A)[[1]] <- list(Aindices)[[1]]
	dimnames(submodel$A)[[2]] <- list(Aindices)[[1]]
	submodel$S <- model$S[Aindices,Aindices]
	dimnames(submodel$S)[[1]] <- list(Aindices)[[1]]
	dimnames(submodel$S)[[2]] <- list(Aindices)[[1]]
	submodel$F <- model$F[selectManifests,Aindices]
	dimnames(submodel$F)[[1]] <- list(selectManifests)[[1]]
	dimnames(submodel$F)[[2]] <- list(Aindices)[[1]]
	if(!is.null(submodel$M)){
		submodel$M <- model$M[1,Aindices]
		submodel$M@values <- t(submodel$M@values)
		submodel$M@labels <- t(submodel$M@labels)
		submodel$M@free <- t(submodel$M@free)
		submodel$M@lbound <- t(submodel$M@lbound)
		submodel$M@ubound <- t(submodel$M@ubound)
	}
	dimnames(submodel$F)[[2]] <-list(Aindices)[[1]]
	if(model@data@type == "raw"){
		submodel$data <- NULL
#		submodel@data@observed <- as.matrix(submodel@data@observed[,selectManifests])
	} else if (model@data@type == "cov") {
		submodel@data@observed <- as.matrix(submodel@data@observed[selectManifests,selectManifests])
#		submodel@data@numObs <- submodel@data@numObs / 2
		colnames(submodel$data@observed) <- selectManifests
		rownames(submodel$data@observed) <- selectManifests
		if (!single.na(submodel@data@means)){
			submodel@data@means <- t(as.matrix(submodel@data@means[1,selectManifests]))
		}
	}
	return(submodel)
}

# the new mxRun function could look like this
# mxRun <- function(model, ..., intervals = FALSE, silent = FALSE, 
		# suppressWarnings = FALSE, unsafe = FALSE,
		# checkpoint = FALSE, useSocket = FALSE, onlyFrontend = FALSE, 
		# useOptimizer = TRUE,usePPML=TRUE){
	# if(!silent) cat("Running", model@name, "\n")
	# frontendStart <- Sys.time()
	# garbageArguments <- list(...)
	# if (length(garbageArguments) > 0) {
		# stop("mxRun does not accept values for the '...' argument")
	# }
	# if(usePPML){
		 # model <- transFormModelDatappml(model)
	# }
	# runHelper(model, frontendStart, intervals,
		# silent, suppressWarnings, unsafe,
		# checkpoint, useSocket, onlyFrontend, useOptimizer)
# }
