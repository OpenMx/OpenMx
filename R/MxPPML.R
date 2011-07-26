imxTransformModelPPML <- function(model) {
	#browser()
	
	###############
	#CHECK SECTION#
	###############
	#checks are ordered from fast to slow
	
	#is the model a RAM model?
	objective <- model$objective
	if(is.null(objective) || !is(objective, "MxRAMObjective")) {
		return(model)
	}
	
	Aname <- objective@A
	Sname <- objective@S
	Fname <- objective@F
	Mname <- objective@M
	Amatrix <- model[[Aname]]
	Smatrix <- model[[Sname]]
	Fmatrix <- model[[Fname]]
	Mmatrix <- NULL
	if (!single.na(Mname)) {
		Mmatrix <- model[[Mname]]
	}
	constraints <- model@constraints
	
	# Matrices must be MxMatrices
	if(!is(Amatrix, "MxMatrix") || !is(Smatrix, "MxMatrix") || !is(Fmatrix, "MxMatrix") || (!single.na(Mname) && !is(Mmatrix, "MxMatrix"))) {
		return(model)
	}
	
	# BASIC: Model must have manifestvars and latentvars and more manifests than latents
	if ( sum(Fmatrix@values) == 0 || sum(Fmatrix@values) == max(dim(Fmatrix@values))
		|| sum(Fmatrix@values) <= (max(dim(Fmatrix@values)) - sum(Fmatrix@values)) ) {
		return(model)
	}
	
	if(!(model@data@type == "raw" || model@data@type == "cov")) {
		return(model)
	}
	
	#are all regression loadings fixed?
	if(any(Amatrix@free)) {
		return(model)	
	}
	
	# If the model is matrix specified, uses numeric indices instead of dimnames
	# Then splits the objective function in the split section
	manifestVars <- NULL
	latentVars <- NULL
	if (length(model@manifestVars) == 0 && length(model@latentVars) == 0) {
		if (length(model@objective@dims) == 1 && single.na(unique(model@objective@dims)))	# no dims anywhere NOTE: is this necessary?
			return(model)
		for (i in 1:max(dim(Fmatrix@values))) {
			if (any(as.logical(Fmatrix@values[,i])))
				manifestVars <- c(manifestVars, as.numeric(i))
			else
				latentVars <- c(latentVars, as.numeric(i))
		}
	}
	else
	{
		manifestVars <- model@manifestVars
		latentVars <- model@latentVars
	}
	
	#only free parameters in the error matrix are on the diagonal?  IN DEVELOPMENT
	#if (any(Smatrix@free[manifestVars, manifestVars] & diag(rep(FALSE,length(manifestVars))))) {
	#	return(model)
	#}
	
	#do all errors have the same label?
	clabels <- array(0,0)
	manHasVar <- matrix(FALSE,length(manifestVars),1)
	rownames(manHasVar) <- manifestVars
	#gather all labels from the direct errors
	for(manifestVar in manifestVars) {
		#if (Smatrix@values[manifestVar,manifestVar] !=0) {
		if (Smatrix@free[manifestVar,manifestVar]) { # NOTE: Does this make sense?
			manHasVar[manifestVar,1] <- TRUE
			clabels <- append(clabels,Smatrix@labels[manifestVar,manifestVar])
		}		
	}
	
	# Classify latents
	fakeLatents <- array(0,0)
	rootLatents <- array(0,0)
	for(latentVar in latentVars){
		# does the latentVar only have one pointer? => fakeLatent
		# This test must include pointers to manifests AND latents
		#  tmp <- which(Amatrix@values[manifestVars,latentVar] != 0) # changed to include latents as well as manifests
		
		# Get nonzero loadings from the asymmetric matrix
		loads <- which(Amatrix@values[,latentVar] != 0)
		# If only 1 loading, might be a fakeLatent
		if (length(loads) == 1) {
			# Can only be a fakeLatent if the one loading is on a manifest var
			if (length(which(Amatrix@values[manifestVars,latentVar] != 0)) > 0) {
				# At this point, latent has to be a fakeLatent or PPML aborts
				# If it points to a manifest and the loading is not 1, or if the manifest
				# already has a variance, the model is not a valid model for PPML
				if ((Amatrix@values[loads,latentVar] == 1) && !manHasVar[loads,1]) {
					fakeLatents <- append(fakeLatents,latentVar)
					next; # Root and Fake categories are mutually exclusive
				} else {
					return(model)
				}
			}
		}
		
		# Check if it's a root latent - no loadings from other latents
		inLoads = which(Amatrix@values[latentVar, latentVars] != 0)
		if (length(inLoads) == 0) {
			rootLatents <- append(rootLatents, latentVar)
			next;
		}
		
		#tmp <- which(Amatrix@values[,latentVar] != 0)	
		#if (length(tmp) == 1) {
		#	#the regression weight has to be one
		#	if ((Amatrix@values[tmp,latentVar] == 1) && !manHasVar[tmp,1]) {
		#		fakeLatents <- append(fakeLatents,latentVar)
		#	} else {
		#		return(model)
		#	}
		#}
	}
	realLatents <- latentVars[is.na(pmatch(latentVars, fakeLatents))]
	nonRootLatents <- latentVars[is.na(pmatch(latentVars, c(fakeLatents, rootLatents)))]
	
	#gather all labels from the fakeLatents
	for(fakeLatent in fakeLatents){
		clabels <- append(clabels, Smatrix@labels[fakeLatent,fakeLatent])
	}
	
	# -------------------------------------------------------------------------
	# Nonhomogenous error:
	# If Cerr is not proportional to I, a transformation is required before applying PPML
	# Cerr is built from the manifest variables in Smatrix and fakeLatents
	
	
	# Rinv will be applied later regardless -- By default, it is an identity matrix
	Rinv <- diag(1, length(manifestVars))
	
	# CHECKING: Any free values off the Cerr diagonal, or the Cerr diagonal is not homogeneous
	if ( any(as.logical(Smatrix@free[manifestVars,manifestVars] & !diag(rep(TRUE, length(manifestVars)))) ) ||
		(length(unique(diag(Smatrix@values[manifestVars, manifestVars]))) > 1) ) {
	
		if (length(fakeLatents) > 0) {
			# Rebuild the model, folding fakeLatents in to the S matrix
			for (fakeLatent in fakeLatents) {
				forWhich <- which(Amatrix@values[,fakeLatent])
				Smatrix[forWhich, forWhich] <- Smatrix@values[fakeLatent,fakeLatent]
			}
			# Remove fakeLatents from A, S, and F matrices
			remainingVars <- c(manifestVars, realLatents)
			Amatrix <- Amatrix[remainingVars, remainingVars]
			Smatrix <- Smatrix[remainingVars, remainingVars]
			Fmatrix <- Fmatrix[ ,remainingVars]
			if (!single.na(Mname)) {
				Mmatrix <- Mmatrix[remainingVars]
			}
			# no fakeLatents remain
			fakeLatents <- array(0,0)
			latentVars <- realLatents
		}
		
		bCERetVal <- buildCErr(Smatrix, manifestVars, constraints)
		Cerr <- bCERetVal[[1]]
		relevantConstraints <- bCERetVal[[2]]
		
		# R <- t(chol(Cerr))  # OLD
		# R <- chol(Cerr)
		# Find transformation matrix R^-1
		#browser()
		Rinv <- solve(chol(Cerr))
		# Apply it to original Cerror to get transformed C'error
		#CerrPrime <- diag(diag(t(Rinv) %*% Cerr %*% Rinv)) # diags to clean up small values
		CerrPrime <- diag(rep(1, dim(Cerr)[1]))
		Smatrix@values[manifestVars, manifestVars] <- CerrPrime # Reinsert to Smatrix
		#Smatrix@free[manifestVars, manifestVars] <- matrix(FALSE, length(manifestVars), length(manifestVars)) | diag(rep(TRUE,length(manifestVars)))
		# Set all labels for the errors to the same label
		Smatrix@free[manifestVars, manifestVars] <- FALSE
		Smatrix@labels[manifestVars, manifestVars] <- NA
		for (manifestVar in manifestVars) {
			if (Smatrix@values[manifestVar, manifestVar] != 0) {
				Smatrix@labels[manifestVar, manifestVar] <- '_PPML_NHEV_ErrParam' #unique(clabels)
				Smatrix@free[manifestVar, manifestVar] <- TRUE # NOTE: ? maybe
			}
		}
		Amatrix[manifestVars, realLatents]@values <- t(Rinv) %*% Amatrix[manifestVars, realLatents]@values # Transform structure matrix
		# Remove relevant constraints
		
		for (relevantConstraint in relevantConstraints) {
			constraints <- setdiff(constraints, list(model@constraints[[relevantConstraint]]))
		}
	}

	###################
	#TRANSFORM SECTION#
	###################
	#get loadings from latents to manifest from A,S,F
	#transform A to E
	
	##browser()
	E <- solve(diag(nrow(Amatrix)) - Amatrix@values)
	#select only these columns of A which are real latents
	#get loadings matrix:= The part of A which goes from the latents to the manifests
	lambda <- as.matrix(E[manifestVars, realLatents])
	qrDecom <- qr(lambda)
	#check if loadings matrix is of full rank
	k <- length(realLatents)
	if (qrDecom$rank != dim(lambda)[2] || k != qrDecom$rank){
		return(model)
	}
	#last check passed, can start modifying model

	#ensure that every free parameter has a unique label
	pair <- omxNameAnonymousParameters(model)
	model <- pair[[1]]
	oldLabels <- pair[[2]]

	#calculate upper triangle matrix
	#orthogonal
	Q <- t(qr.Q(qrDecom, complete = TRUE))
	#transform loadings matrix
	lambda <- Q %*% lambda
	#set all rows from k+1 to zero
	lambda[(k + 1) : nrow(lambda), ] <- 0
	#insert new loadings matrix in A
	Amatrix@values[manifestVars,realLatents] <- lambda
	Amatrix@values[realLatents,realLatents] <- 0 # For latents predicting latents
	model[[Aname]] <- Amatrix
	#insert modified S, F, M matrices (from Matrix Spec, NHEV) in to model
	model[[Sname]] <- Smatrix
	model[[Fname]] <- Fmatrix
	if (!single.na(Mname)) {
		model[[Mname]] <- Mmatrix
	}
	if (!is.numeric(latentVars)) {
		model@latentVars <- latentVars
	} else if (length(latentVars) < (max(dim(Fmatrix@values)) - sum(Fmatrix@values)) ) {
		# if number of latentVars has been reduced, remove appropriate variable
		# names from the objective
		oldLatentVars <- NULL
		for (i in 1:max(dim(Fmatrix@values))) {
			if (!any(as.logical(Fmatrix@values[,i])))
				oldLatentVars <- c(latentVars, as.numeric(i))
		}
		model@objective@dims <- model@objective@dims[-setdiff(oldLatentVars, latentVars)]
	}
	model@constraints <- constraints
	
	
	
	#tranform data or cov
	if(model$data@type == "raw") {
		model$data@observed <- as.matrix(model@data@observed) %*% (Rinv) %*% t(Q)
		colnames(model$data@observed) <- manifestVars
	}
	else if(model$data@type == "cov") {
		# Transform variance data
		model$data@observed <- (Q) %*% t(Rinv) %*% as.matrix(model@data@observed) %*% (Rinv) %*% t(Q)
		
		# Restore dimnames
		if (!(is.numeric(manifestVars) && is.numeric(latentVars))) {
			# For path-specified models:
			colnames(model$data@observed) <- manifestVars
			rownames(model$data@observed) <- manifestVars
		} else {
			# For matrix-specified models
			colnames(model$data@observed) <- model@objective@dims[manifestVars]
			rownames(model$data@observed) <- model@objective@dims[manifestVars]
		}
		
		# Transform means data, if it exists
		if (!single.na(model@data@means)){
			model$data@means <- t(Q %*% t(as.matrix(model@data@means)))
			# Restore dimnames
			if (!(is.numeric(manifestVars) && is.numeric(latentVars)))
				colnames(model$data@means) <- manifestVars	# Path-specified
			else
				colnames(model$data@means) <- model@objective@dims[manifestVars] # Matrix-specified
		}
		
	}
	
	###############
	#SPLIT SECTION#
	###############
	#leftmodel
	#calculate A realLatents indices for left model
	#first k manifest vars
	selectManifests <- manifestVars[1:k]
	#all "REAL" latents + fake latents which point into the selectManifests
	#remember that we checked that each fakeLatent column in A only has 1 non zero entry
	#which is one
	indices <- colSums(Amatrix@values[selectManifests, fakeLatents, drop = FALSE])

	selectFake <- fakeLatents[which(indices == 1)]
	selectLatents <- append(realLatents, selectFake)	
	leftmodel <- mxRename(model, 'leftmodel')	
	leftmodel <- selectSubModelFData(leftmodel, selectLatents, selectManifests)

	
	#rightmodel
	selectLatents <- latentVars[is.na(pmatch(latentVars, selectLatents))]
	selectManifests <- manifestVars[is.na(pmatch(manifestVars, selectManifests))]
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

#@author:	Julian Karch
#		jk3nq@virginia.edu
#@idea>		Timo von Oertzen
#		timo@virginia.edu
selectSubModelFData <- function(model, selectLatents, selectManifests) {
	Aindices <- append(selectManifests, selectLatents)
	
	#build the new model
	submodel <- model
	
	submodel$A <- model$A[Aindices,Aindices]
	submodel$S <- model$S[Aindices,Aindices]
	submodel$F <- model$F[selectManifests,Aindices]
	
	# Restore dimnames to new matrices
	# Only necessary for path specified models
	if (!(is.numeric(selectManifests) && is.numeric(selectLatents))) {
		submodel@manifestVars <- selectManifests
		submodel@latentVars <- selectLatents
		dimnames(submodel$A)[[1]] <- list(Aindices)[[1]]
		dimnames(submodel$A)[[2]] <- list(Aindices)[[1]]
		dimnames(submodel$S)[[1]] <- list(Aindices)[[1]]
		dimnames(submodel$S)[[2]] <- list(Aindices)[[1]]
		dimnames(submodel$F)[[1]] <- list(selectManifests)[[1]]
		dimnames(submodel$F)[[2]] <- list(Aindices)[[1]]
	}
	else
	{
		submodel@objective@dims <- submodel@objective@dims[Aindices]
	}
	
	if(!is.null(submodel$M)){
		submodel$M <- model$M[1,Aindices]
		# submodel$M@values <- t(submodel$M@values)
		# submodel$M@labels <- t(submodel$M@labels)
		# submodel$M@free <- t(submodel$M@free)
		# submodel$M@lbound <- t(submodel$M@lbound)
		# submodel$M@ubound <- t(submodel$M@ubound)
	}
	
	if(model@data@type == "raw"){
		submodel$data <- NULL
#		submodel@data@observed <- as.matrix(submodel@data@observed[,selectManifests])
	} else if (model@data@type == "cov") {
		# Pull out the proper manifest variables from the covariance matrix
		# to create the appropriate covariance matrix for the submodel
		submodel@data@observed <- as.matrix(submodel@data@observed[selectManifests,selectManifests])
	
		# Restore dimnames to the new covariance matrix, if path specified
		if (!(is.numeric(selectManifests) && is.numeric(selectLatents))) {
			colnames(submodel$data@observed) <- selectManifests
			rownames(submodel$data@observed) <- selectManifests
		}
		# If means data exists, pull out the appropriate manifest variables
		if (!single.na(submodel@data@means)){
			submodel@data@means <- t(as.matrix(submodel@data@means[1,selectManifests]))
		}
	}
	return(submodel)
}

buildCErr <- function(Smatrix, manifestVars, constraints) {
	# Build Cerr
	# Get labels
	errLabels <- setdiff(unique(as.vector(Smatrix@labels[manifestVars,manifestVars])), NA)
	Cerr <- matrix(data=0, nrow=length(manifestVars), ncol=length(manifestVars), dimnames=list(manifestVars, manifestVars))
	
	# Label[1] * Factor == Label[2]
	conLabels <- as.character(array(0,0))
	conFactors <- as.numeric(array(0,0))
	
	relevantConstraints <- as.character(array(0,0))
	for (constraint in constraints) {
		# NOTE: Necessary to get a length check on the constraint?
		
		# Constraint is relevant if left and right hand side each involve
		# one label from the Smatrix, and no other labels
		# If a constraint is a relationship between terms in the Smatrix
		# and terms outside of the Smatrix, then the model is not of the
		# correct form.
		# NOTE: This would not be true in the case of a label on a fixed value
		
		hasRelevantLabelLeft <- FALSE
		hasRelevantLabelRight <- FALSE
		hasIrrelevantLabelLeft <- FALSE
		hasIrrelevantLabelRight <- FALSE
		# Left side
		# NOTE: Constraint checking algorithm should only return(model) if the
		# constraint involves at least one relevant label. Otherwise, it can 
		# just ignore the constraint.
		for (term in 1:length(constraint@formula[[2]])) {
			# If the term is a label, do checks
			if (is.character(constraint@formula[[2]][[term]]) ) {
				if (any(as.character(constraint@formula[[2]][[term]]) == errLabels)) {
					if (hasRelevantLabelLeft || hasIrrelevantLabelLeft) {
						# Only one label per side allowed
						return(model)
					}
					hasRelevantLabelLeft <- TRUE
				} else {
					hasIrrelevantLabelLeft <- TRUE
					if (hasRelevantLabelLeft) {
						# Can't have irrelevant and relevant labels in the same constraint
						return(model)
					}
				}
			}
		}
		# Right side
		for (term in 1:length(constraint@formula[[3]])) {
			# If the term is a label, do checks
			if (is.character(constraint@formula[[3]][[term]]) ){
				if (any(as.character(constraint@formula[[3]][[term]]) == errLabels)) {
					if (hasIrrelevantLabelLeft) {
						# Can't have an irrelevant label on left and relevant on right
						return(model)
					}
					if (hasIrrelevantLabelRight || hasRelevantLabelRight) {
						# Can't have more than one label on a side
						return(model)
					}
					hasRelevantLabelRight <- TRUE
				} else {
					hasIrrelevantLabelRight <- TRUE
					if (hasRelevantLabelLeft) {
						# Can't have a relevant label left and an irrelevant right
						return(model)
					}
					if (hasRelevantLabelRight) {
						# Can't have two labels on one side
						return(model)
					}
				}
			}
		}
		
		# If not relevant, continue to next constraint
		if (!(hasRelevantLabelLeft && hasRelevantLabelRight)) {
			next
		} else {
			relevantConstraints <- c(relevantConstraints, constraint@name)
		}
		
		# If the constraint is not an equality constraint, the model
		# is not in the correct form
		if (constraint@formula[[1]] != '==') {
			return(model)
		}
		
		# Left and right hand sides must be of length 3 or 1
		# NOTE: Not a very valuable test, necessary at all?
		# NOTE: Could fold in to below loop
		if (!( (length(constraint@formula[[2]]) == 1 || length(constraint@formula[[2]]) == 3) &&
			(length(constraint@formula[[3]]) == 1 || length(constraint@formula[[3]]) == 3) ) ) {
			return(model)
		}
		
		newFactor <- 1
		newLabel <- c("", "")
		# examine left, then right side
		for (side in 2:3) {
			if (length(constraint@formula[[side]]) == 1) {
				# Earlier check ensured if there is only one term on a side,
				# it's a relevant label
				newLabel[side-1] <- as.character(constraint@formula[[side]])
			} else {
				# Length 3
				if (constraint@formula[[side]][[1]] != '*') {
					# Relationship must be multiplicative
					return(model)
				}
				if (is.character(constraint@formula[[side]][[2]]) && is.numeric(constraint@formula[[side]][[3]]) ){
					newLabel[side-1] <- as.character(constraint@formula[[side]][[2]])
					if (side == 1) {
						newFactor <- newFactor * as.numeric(constraint@formula[[side]][[3]])
					} else {
						newFactor <- newFactor / as.numeric(constraint@formula[[side]][[3]])
					}	
				} else if (is.character(constraint@formula[[side]][[3]]) && is.numeric(constraint@formula[[side]][[2]]) ){
					newLabel[side-1] <- as.character(constraint@formula[[side]][[3]])
					if (side == 1) {
						newFactor <- newFactor * as.numeric(constraint@formula[[side]][[2]])
					} else {
						newFactor <- newFactor / as.numeric(constraint@formula[[side]][[2]])
					}
				} else {
					# Something is wrong, not a numeric times a label
					# NOTE: this will occur in cases where you have the label multiplied by
					# the product of two numeric factors (i.e., 0.1*0.2*Err)
					return(model)
				}
			}
		}
		conLabels <- rbind(conLabels, newLabel)
		conFactors <- c(conFactors, newFactor)
	}

	errValues <- rep(0, length(errLabels))
	errValueSet <- rep(FALSE, length(errLabels))
	errValues[1] <- 1
	errValueSet[1] <- TRUE
	errConUsed <- rep(FALSE, dim(conLabels)[1])
	
	for (i in 1:length(errLabels)) {
		if (errValueSet[i]) {
			# propagate out, looking for inconsistencies
			for (j in 1:dim(conLabels)[1]) {
				if (!errConUsed[j]) {
					if (conLabels[j,1] == errLabels[i]) {
						otherLabel <- which(errLabels == conLabels[j,2])
						if (errValueSet[otherLabel]) {
							if (errValues[otherLabel] != errValues[i] * conFactors[j]) {
								# inconsistency
								return(model)
							}
							errConUsed[j] <- TRUE
						} else {
							errValues[otherLabel] <- errValues[i] * conFactors[j]
							errValueSet[otherLabel] <- TRUE
							errConUsed[j] <- TRUE
						}
					} else if (conLabels[j,2] == errLabels[i]) {
						otherLabel <- which(errLabels == conLabels[j,1])
						if (errValueSet[otherLabel]) {
							if (errValues[otherLabel] != errValues[i] / conFactors[j]) {
								# inconsistency
								return(model)
							}
							errConUsed[j] <- TRUE
						} else {
							errValues[otherLabel] <- errValues[i] / conFactors[j]
							errValueSet[otherLabel] <- TRUE
							errConUsed[j] <- TRUE
						}
					}
				}
			}
		} else {
			# propagate in, looking for inconsistences
			for (j in 1:dim(conLabels)[1]) {
				if (!errConUsed[j]) {
					if (conLabels[j,1] == errLabels[i]) {
						otherLabel <- which(errLabels == conLabels[j,2])
						if (errValueSet[otherLabel]) {
							if (errValueSet[i] && (errValues[i] != (errValues[otherLabel] / conFactors[j])) ) {
								# inconsistency detection after initially set
								return(model)
							}
							errValues[i] <- errValues[otherLabel] / conFactors[j]
							errValueSet[i] <- TRUE
							errConUsed[j] <- TRUE
						}
					} else if (conLabels[j,2] == errLabels[i]) {
						otherLabel <- which(errLabels == conLabels[j,1])
						if (errValueSet[otherLabel]) {
							if (errValueSet[i] && (errValues[i] != (errValues[otherLabel] * conFactors[j])) ) {
								# inconsistency detection after initially set
								return(model)
							}
							errValues[i] <- errValues[otherLabel] * conFactors[j]
							errValueSet[i] <- TRUE
							errConUsed[j] <- TRUE
						}
					}
				}
			}
		}
	}
	
	if (any(!errValueSet)) {
		# Unconstrained free parameter in the symmetric matrix, model is not of the
		# proper form
		return(model)
	}
	
	# Create Cerr
	for (errLabel in errLabels) {
		Cerr[Smatrix@labels[manifestVars, manifestVars] == errLabel] <- errValues[which(errLabels == errLabel)]
	}
	
	return(list(Cerr, relevantConstraints))
}

imxRestoreResultPPML <- function(model, result) {
	paramLabels <- names(result@output$estimate)		# Get list of parameter labels (some are unknown***, NA)
	paramValues <- as.vector(result@output$estimate)	# Get corresponding list of values
	
	### Set parameters working around omxSetParameters not liking unnamed free params
	# Free parameters are always in the S matrix and sometimes in the M matrix
	# Values with labels can be in multiple positions, set them first - omxSetParameters can handle this
	
	# Check if a transform from the NHEV case has been performed
	# If so, find the error covariance matrix and find the value of each
	# parameter from the error scaling parameter and the error matrix
	if (!is.na(any(paramLabels == "_PPML_NHEV_ErrParam")) && any(paramLabels == "_PPML_NHEV_ErrParam")) {
		# pull out Smatrix
		Smatrix <- model@matrices$S
		
		# Check for fakeLatents
		fakeLatents <- list()
		for (latentVar in model@latentVars) {
			if (sum(model@matrices$A@values[,latentVar]) == 1) {
				if (sum(model@matrices$A@values[model@manifestVars, latentVar]) == 1) {
					fakeLatents <- c(fakeLatents, latentVar)
				}
			}
		}
		
		if (length(fakeLatents) > 0) {
			# Rebuild the model, folding fakeLatents in to the S matrix
			Smatrix@values[manifestVars, manifestVars] <- Cerr # Insert Cerr in to Smatrix
			# Remove fakeLatents from A, S, and F matrices
			remainingVars <- c(manifestVars, realLatents)
			Amatrix <- Amatrix[remainingVars, remainingVars]
			Smatrix <- Smatrix[remainingVars, remainingVars]
			Fmatrix <- Fmatrix[,remainingVars]
			# no fakeLatents remain
			fakeLatents <- array(0,0)
			latentVars <- realLatents
		}		
		
		Cerr <- buildCErr(Smatrix, model@manifestVars, model@constraints)
		
		errScale <- paramValues[which(paramLabels == "_PPML_NHEV_ErrParam")]
		
		newParamLabels <- array(0,0)
		newParamValues <- array(0,0)
		for (errParam in setdiff(unique(model@matrices$S@labels), c(paramLabels, NA))) {
			newParamLabels <- c(newParamLabels, errParam)
			newParamValues <- c(newParamValues, errScale * setdiff(unique(Cerr[[1]][which(Smatrix@labels[model@manifestVars, model@manifestVars] == errParam)]), NA))
		}
		PNEPindex <- which(paramLabels == "_PPML_NHEV_ErrParam")
		paramValues <- c(paramValues[-PNEPindex], newParamValues)
		paramLabels <- c(paramLabels[-PNEPindex], newParamLabels)
	}
	
	allSLabels = unique(model@matrices$S@labels[which(!is.na(as.vector(model@matrices$S@labels)))])	# Get list of labels in S
	
	labeledSValues = NULL;
	labeledSLabels = NULL;
	if (length(allSLabels) > 0) {
		# Find location of all the S labels in paramLabels
		sLabelLIV <- rep(FALSE, length(paramLabels)) # First build a logical index vector
		for (sLabel in allSLabels) {
			sLabelLIV <- sLabelLIV | (paramLabels == sLabel)
		}
		sLabelLoc <- which(sLabelLIV)	# Then get the indices
				
		labeledSValues <- paramValues[sLabelLoc]	# Extract values of labeled params
		labeledSLabels <- paramLabels[sLabelLoc]	# Create corresponding label list
	}
	
	# Get labeledMValues, labeledMLabels if the M matrix exists
	allMLabels = NULL;
	labeledMValues = NULL;
	labeledMLabels = NULL;
	if (!is.null(model@matrices$M)) {
		allMLabels = unique(model@matrices$M@labels[which(!is.na(as.vector(model@matrices$M@labels)))])	# Get list of labels in M	
		
		if (length(allMLabels) > 0) {
			# Find location of all the M labels in paramLabels
			mLabelLIV <- rep(FALSE, length(paramLabels)) # First build a logical index vector
			for (mLabel in allMLabels) {
				mLabelLIV <- mLabelLIV | (paramLabels == mLabel)
			}
			mLabelLoc <- which(mLabelLIV) # Then get the indices
			
			labeledMValues <- paramValues[mLabelLoc]	# Extract values of labeled params
			labeledMLabels <- paramLabels[mLabelLoc] # Create corresponding label list
		}
	}

	# Combine all the labeled values together for omxSetParameters to use
	labeledValues = c(labeledSValues, labeledMValues)
	labeledLabels = c(labeledSLabels, labeledMLabels)
	
	if (length(labeledValues) > 0) {	# Make sure there are labeled values
		model <- omxSetParameters(model, labels = labeledLabels, values = labeledValues, strict = FALSE)	# Set labeled params
	}
	
	#browser()
	
	## Unlabeled parameters need to be set manually within the matrices
	# They are always in the order they appear in the S matrix, and then the M matrix
	allLabels <- c(allSLabels, allMLabels)
	if ( (length(paramValues) - length(allLabels)) > 0) {	# if there are any unlabeled params
	
		# Find location of all the unlabeled params in paramLabels
		unlabeledLIV <- rep(TRUE, length(paramLabels)) # First build a logical index vector
		for (label in allLabels) {
			unlabeledLIV <- unlabeledLIV & (paramLabels != label)
		}
		unlabeledLIV[which(is.na(unlabeledLIV))] <- TRUE
		unlabeledLoc <- which(unlabeledLIV) # Then get the indices
		
		unlabeledValues <- paramValues[unlabeledLoc]	# Extract unlabeled values
		
		# There must be the correct number of free parameter values in unlabeledValues or something has gone horribly wrong
		# Set unlabeled S values - Must set these first as S values are first in order in unlabeledValues
		# In order to deal with off-diagonal values, set values on the upper triangular matrix, and then rebuild symmetry
		paramIndicesS = which( as.vector(model@matrices$S@free & upper.tri(model@matrices$S@free, diag=TRUE)) 
			& as.vector(is.na(model@matrices$S@labels)) )	# Find indices of free unlabeled S values
		if (length(paramIndicesS) > 0) {	# Make sure there are unlabeled parameters in S
			model@matrices$S@values[paramIndicesS] <- unlabeledValues[1:length(paramIndicesS)]		# Plug in unlabeled S values
			# And restore symmetry
			model@matrices$S@values[which(lower.tri(model@matrices$S@values))] <- t(model@matrices$S@values)[which(lower.tri(model@matrices$S@values))]
		}
		
		# If there is an M matrix, plug the param values in to that
		if (!is.null(model@matrices$M)) {
			paramIndicesM = which( as.vector(model@matrices$M@free) & as.vector(is.na(model@matrices$M@labels)) )	# Find indices of free unlabeled M values
			if (length(paramIndicesM) > 0) {	# Make sure there are unlabeled parameters in M
				model@matrices$M@values[paramIndicesM] <- unlabeledValues[(length(paramIndicesS)+1):length(unlabeledValues)]	# Plug in unlabeled M values
			}
		}
	}
	
	# Generate fixed model by running with the parameters set to the correct values on the 
	# original model without using the optimizer
	fixed <- mxRun(model, useOptimizer = FALSE, silent = TRUE)
	
	# Transfer important data from original result:
	# Add elapsed times together
	fixed@output$frontendTime <- fixed@output$frontendTime + result@output$frontendTime
	fixed@output$backendTime <- fixed@output$backendTime + result@output$backendTime
	fixed@output$independentendTime <- fixed@output$independentendTime + result@output$independentendTime
	fixed@output$wallTime <- fixed@output$wallTime + result@output$wallTime
	fixed@output$cpuTime <- fixed@output$cpuTime + result@output$cpuTime
	
	# Carry over statuses
	fixed@output$status <- result@output$status
	
	# Model is un-PPMLed and relevant data has been transferred, return
	return(fixed)	
}

single.na <- function(a) {
    return((length(a) == 1) &&
        (is.list(a) || is.vector(a) || is.matrix(a)) &&
        (is.na(a) == TRUE))
}	

# DEVELOPMENT FUNCTIONS

# Function could be extended to be more general, but for now is just for comparing fits for models,
# their PPML transforms, and the corresponding reversed PPML transforms
imxCheckFitsPPML <- function(res1, res2, checkHessians = TRUE, checkLL = TRUE) {
	# Check -2logLLs versus each other
	if (checkLL)
		omxCheckCloseEnough(res1@output$Minus2LogLikelihood, res2@output$Minus2LogLikelihood, 0.001)
	
	# Check parameters versus each other
	# Parameters will be in the same order, so a simple iteration should work
	for ( i in 1:length(res1@output$estimate) ) {
		omxCheckCloseEnough(res1@output$estimate[i], res2@output$estimate[i], 0.001)
	}

	# Hessians will be wrong for reversed PPML transforms
	if (checkHessians) {
		# Very loose epsilons for these
		# Check hessianCholeskys versus each other
		omxCheckCloseEnough(res1@output$hessianCholesky, res2@output$hessianCholesky, 0.3)
		
		# Check calculatedHessians versus each other
		omxCheckCloseEnough(res1@output$calculatedHessian, res2@output$calculatedHessian, 0.3)	
		
		# Estimated Hessians are way, way far off from each other -- doesn't seem worth checking
		## Check estimatedHessians versus each other
		##omxCheckCloseEnough(res1@output$estimatedHessian, res2@output$estimatedHessian, 0.5)	
	}
}

jlRun <- function(model){
	model <- imxTransformModelPPML(model)
	return(mxRun(model))
}

dhRun <- function(model) {
	PPMLmodel <- imxTransformModelPPML(model)
	result <- mxRun(PPMLmodel)
	return(imxRestoreResultPPML(model, result))
}

imxTestPPML <- function(model, checkLL = TRUE) {
	#browser()
	res1 <- mxRun(model) # Standard fit
	res2 <- mxRun(imxTransformModelPPML(model))	# PPML fit
	res3 <- imxRestoreResultPPML(model, res2)	# Reverse transform PPML
	#browser()
	
	# NOTE: Not checking Hessians anymore, they're never very close -- not sure that
	# a comparison is actually meaningful
	# NOTE: checkLL parameter can turn off log likelihood checking for this test,
	# useful for non-homogeneous error variance cases where the transformed log 
	# likelihood is different.
	imxCheckFitsPPML(res1, res2, checkLL = checkLL, checkHessians = FALSE) # Check standard fit vs PPML model
	
	# NOTE Always check likelihood for this test, to make sure everything was plugged in correctly
	imxCheckFitsPPML(res1, res3, checkLL = TRUE, checkHessians = FALSE)	# Check standard fit vs reverse transformed PPML model
	
	# NOTE This last check might be mostly unnecessary.  mxCheckFitPPML with checkHessians = FALSE currently
	# only checks the parameters and -2logLL. The parameters are copied over from res2; thus, this
	# check effectively only compares the -2logLLs between res1 and res3 and then wastes time rechecking
	# the parameters.
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
