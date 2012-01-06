imxCheckPPMLApplicable <- function(model) {
	# browser()
	###############
	#CHECK SECTION#
	###############
	#checks are ordered from fast to slow
	
	solveType <- "Solve"
	
	#is the model a RAM model?
	# NOTE: This could be made redundant by implementing the transform call
	# from the MxRAMObjective function
	# Depends on FIML implementation
	objective <- model$objective
	if(is.null(objective) || !is(objective, "MxRAMObjective")) {
		return(NA)
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
		return(NA)
	}
	
	if( !is.null(model@data) ) {
		if ( !(model@data@type == "raw" || model@data@type == "cov") )  {
			return(NA)
		}
	} else {
		# TODO: If the model@data is NULL, optimized solution is still possible!
		#solveType = "Split"
		return(NA)
	}
	
	#are all regression loadings fixed?
	if(any(Amatrix@free)) {
		return(NA)	
	}
	
	# If the model is matrix specified, uses numeric indices instead of dimnames
	# Then splits the objective function in the split section
	manifestVars <- NULL
	latentVars <- NULL
	if (length(model@manifestVars) == 0 && length(model@latentVars) == 0) {
		if (length(model@objective@dims) == 1 && single.na(unique(model@objective@dims)))	# no dims anywhere NOTE: is this necessary?
			return(NA)
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
		
	# Classify latents
	fakeLatents <- array(0,0)
	rootLatents <- array(0,0)
	for(latentVar in latentVars){
		# Get nonzero loadings from the asymmetric matrix
		loads <- which(Amatrix@values[,latentVar] != 0)
		
		# All loadings from the latent must be on manifestVars
		#  TODO: Latents that covary with other latents but have no loadings on to anything?
		#  Lambda is probably singular in these cases.
		if (length(which(Amatrix@values[latentVars, latentVar] != 0)) == 0) { # No loadings on latentVars
			# Latent must not covary with any manifests or other latents
			if ( Smatrix@free[latentVar,latentVar] && length(which(Smatrix@free[,latentVar])) == 1) {
				# All of these manifestVars must have no inherent error variances /
				# loads are all on to manifests without their own error variances
				#  TODO: Overlapping fakelatents; keep track of which manifests' error variances have
				#  already been "claimed" by a fakelatent.  This case gets complicated because we would
				#  have to make a choice about which latent to call the error variance and which to leave
				#  as a latent.  Criteria for best choice: consistency with error variance labels,
				#  consistency with NHEV constraints, etc.
				if ( all(!diag(Smatrix@free[loads, loads])) && all(Smatrix@values[loads, loads] ==0) ) {
					# Loadings must all have the same value
					if ( length(unique(Amatrix@values[loads, latentVar])) == 1 ) {
						# DEVELOPMENT: Loadings must all have value 1
						#  TODO: This is not necessary as long as they're all the same,
						#  just makes folding the fakelatents back out harder!
						if ( all(Amatrix@values[loads, latentVar] == 1) ) {
							fakeLatents <- append(fakeLatents,latentVar)
							next; # Root and Fake categories are mutually exclusive (explicitly, not by property; this wouldn't work without this code.)
						}
					}
				}
			}
		}
		
		# Check if it's a root latent - no loadings from other latents
		inLoads = which(Amatrix@values[latentVar, latentVars] != 0)
		if (length(inLoads) == 0) {
			rootLatents <- append(rootLatents, latentVar)
		}
	}
	realLatents <- latentVars[is.na(pmatch(latentVars, fakeLatents))]
	nonRootLatents <- latentVars[is.na(pmatch(latentVars, c(fakeLatents, rootLatents)))]
	
	# DEVeLOPMENT -- Some of these models should be transformable
	if (length(nonRootLatents) > 0) {
		return(NA)
	}
	
	# DEVELOPMENT -- Can't handle fakeLatents yet
	if ( length(fakeLatents) > 0 ) {
		return(NA)
	}
	
	# Structure matrix must be invertible
	if ( det( diag(length(realLatents)) - Amatrix@values[realLatents, realLatents] ) == 0 ) {
		return(NA)
	}
	
	# BASIC: Model must have manifestvars and latentvars and more manifests than latents
	if (length(realLatents) == 0 || length(manifestVars) == 0 || 
		length(realLatents) >= length(manifestVars)) {
		return(NA)
	}
	
	# TODO: NHEV CHECK
	# For now, filter out all possible NHEV cases -- should be analytically solvable eventually
	# TODO: In NHEV cases, constraints that aren't part of the error covariance matrix -- can PPML
	# be applied with any other constraint on the model, or will further constraints always break it?
	
	# Fold in fakelatents
	if (length(fakeLatents) > 0) {
		# Rebuild the model, folding fakeLatents in to the S matrix
		for (fakeLatent in fakeLatents) {
			forWhich <- which(Amatrix@values[ ,fakeLatent] != 0) # Manifests to which the fakeLatent corresponds
			if ( is.na(Smatrix@labels[fakeLatent, fakeLatent]) ) {
				# Fake latent must have a non-NA label so that the error variances are properly constrained
				Smatrix@labels[fakeLatent, fakeLatent] <- paste("PPML_FL_", which(fakeLatents == fakeLatent), sep="")
			}
			if (length(forWhich) == 1) {
				Smatrix@values[forWhich, forWhich] <- Smatrix@values[fakeLatent,fakeLatent] # Move starting values
				Smatrix@free[forWhich, forWhich] <- TRUE # The manifests are now free
				Smatrix@labels[forWhich, forWhich] <- Smatrix@labels[fakeLatent,fakeLatent] # Give these manifests the fakeLatent's label				
			} else if (length(forWhich) > 1) {
				diag(Smatrix@values[forWhich, forWhich]) <- rep(Smatrix@values[fakeLatent,fakeLatent], length(forWhich)) # Move starting values
				diag(Smatrix@free[forWhich, forWhich]) <- rep(TRUE, length(forWhich)) # The manifests are now free
				# Give these manifests the fakeLatent's label
				diag(Smatrix@labels[forWhich, forWhich]) <- rep(Smatrix@labels[fakeLatent,fakeLatent], length(forWhich))
			}
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
		# fakeLatents <- array(0,0)
		latentVars <- realLatents
	}
	
	# Any free values off the diagonal in the error matrix?
	if ( any(as.logical(Smatrix@free[manifestVars,manifestVars] & !diag(rep(TRUE, length(manifestVars)))) ) ) {
		# DEVELOPMENT: For now, reject any NHEV case
		return(NA)
	}	
	# Any constraints/multiple diagonal labels/no diagonal labels -> NHEV/NA
	if (length(model@constraints) > 0 || length(unique(diag(Smatrix[manifestVars,manifestVars]@labels))) > 1
		|| all(is.na(Smatrix[manifestVars,manifestVars]@labels)) ) {
		# DEVELOPMENT: For now, reject any NHEV case
		return(NA)
	}
	
	# DATA MISSINGNESS CHECK	
	if (model@data@type == "raw" && any(is.na(model@data@observed))) {
		# Check missingness patterns: if none of the patterns have more manifests than latents,
		# then PPML cannot be applied at all
		if ( (length(manifestVars) - min(apply(is.na(model$data@observed),1,sum))) <= length(realLatents) ) {
			return(NA)
		}
		# Check missingness patterns: if any of them have fewer manifests than the model has latents,
		# then the solution cannot be done analytically
		# TODO:  If there are patterns where a latent has a zero loading on to all of the manifests,
		# can remove that latent and that pattern might be PPML-transformable after all
		if ( (length(manifestVars) - max(apply(is.na(model$data@observed),1,sum))) >= length(realLatents) ) {
			solveType <- "Split"
		}
		# solveType <- "Solve"
		solveType <- "Split" # TODO: SHOULD be analytical, but not functioning yet
	}
		
	# Partial Solve possible when all covariances between latents are fixed to zero but all
	# variances of the latents are free
	# All covariances are fixed, zero
	if ( !any(Smatrix@free[realLatents, realLatents] & !diag(TRUE, length(realLatents), length(realLatents)) )
		&& all( abs(Smatrix@values[realLatents, realLatents] - diag(diag(Smatrix@values[realLatents,realLatents]))) < .001 ) ) {
		if ( all(diag(Smatrix@free[realLatents, realLatents])) ) {
			if (solveType != "Split") solveType <- "PartialSolve"
		} else {
			solveType <- "Split"
		}
	} else if ( !all(Smatrix@free[realLatents, realLatents]) ) {
		return(NA)
	}
	# Although not in the code,
	# solveType <- "Solve" happens implicitly when Smatrix@free[latentVars, latentVars] is all free
	# (saturated covariance matrix)
	
	return(solveType)
	
}

imxTransformModelPPML <- function(model, solveType = "Check") {
	
	# browser()
	
	### CHECK SECTION
	### Uses imxCheckPPMLApplicable
	
	# Check if a non-default value hasn't been provided
	# NOTE: Allowing a value to be provided allows the check to be cut out for
	# time savings in cases where a number of models (or the same model with different
	# datasets) with known PPML applicability, are being run in batches.
	if (solveType == "Check") {
		solveType <- imxCheckPPMLApplicable(model)
	}
	
	if (is.na(solveType)) return(model)
	
	# IN-DEVELOPMENT : Only analytical solution implemented currently
	# So, for split and partially solved models, abort the transform
	if (solveType == "PartialSolve" || solveType == "Split") {
	# if (solveType == "Split") { # DEVELOPMENT
		return(model)
	}
	
	# -------------------------------------------------------------------------
	### SET-UP Section
	### Set up for the transform
	objective <- model$objective
	
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
	
	
	# browser()
	# Classify latents
	fakeLatents <- array(0,0)
	rootLatents <- array(0,0)
	for(latentVar in latentVars){
		# Get nonzero loadings from the asymmetric matrix
		loads <- which(Amatrix@values[,latentVar] != 0)
		
		# All loadings from the latent must be on manifestVars
		#  TODO: Latents that covary with other latents but have no loadings on to anything?
		#  Lambda is probably singular in these cases.
		if (length(which(Amatrix@values[latentVars, latentVar] != 0)) == 0) { # No loadings on latentVars
			# Latent must not covary with any manifests or other latents
			if ( Smatrix@free[latentVar,latentVar] && length(which(Smatrix@free[,latentVar])) == 1) {
				# All of these manifestVars must have no inherent error variances /
				# loads are all on to manifests without their own error variances
				#  TODO: Overlapping fakelatents; keep track of which manifests' error variances have
				#  already been "claimed" by a fakelatent.  This case gets complicated because we would
				#  have to make a choice about which latent to call the error variance and which to leave
				#  as a latent.  Criteria for best choice: consistency with error variance labels,
				#  consistency with NHEV constraints, etc.
				if ( all(!diag(Smatrix@free[loads, loads])) && all(Smatrix@values[loads, loads] ==0) ) {
					# Loadings must all have the same value
					if ( length(unique(Amatrix@values[loads, latentVar])) == 1 ) {
						# DEVELOPMENT: Loadings must all have value 1
						#  TODO: This is not necessary as long as they're all the same,
						#  just makes folding the fakelatents back out harder!
						if ( all(Amatrix@values[loads, latentVar] == 1) ) {
							fakeLatents <- append(fakeLatents,latentVar)
							next; # Root and Fake categories are mutually exclusive (explicitly, not by property; this wouldn't work without this code.)
						}
					}
				}
			}
		}
		
		# Check if it's a root latent - no loadings from other latents
		inLoads = which(Amatrix@values[latentVar, latentVars] != 0)
		if (length(inLoads) == 0) {
			rootLatents <- append(rootLatents, latentVar)
		}
	}
	realLatents <- latentVars[is.na(pmatch(latentVars, fakeLatents))]
	nonRootLatents <- latentVars[is.na(pmatch(latentVars, c(fakeLatents, rootLatents)))]
	
	# -------------------------------------------------------------------------
	# IN DEVELOPMENT
	# Fake Latents
	# There are multiple ways to specify the error variance terms.  There is the
	# usual, direct way of allowing the term in the S matrix to be free, but it
	# can also be specified using latent variables.
	# 
	# This segment adjusts the model so that all error variance is specified using
	# only the S matrix, without any latent variables
	if (length(fakeLatents) > 0) {
		# Rebuild the model, folding fakeLatents in to the S matrix
		for (fakeLatent in fakeLatents) {
			forWhich <- which(Amatrix@values[ ,fakeLatent] != 0) # Manifests to which the fakeLatent corresponds
			if ( is.na(Smatrix@labels[fakeLatent, fakeLatent]) ) {
				# Fake latent must have a non-NA label so that the error variances are properly constrained
				Smatrix@labels[fakeLatent, fakeLatent] <- paste("PPML_FL_", which(fakeLatents == fakeLatent), sep="")
			}
			if (length(forWhich) == 1) {
				Smatrix@values[forWhich, forWhich] <- Smatrix@values[fakeLatent,fakeLatent] # Move starting values
				Smatrix@free[forWhich, forWhich] <- TRUE # The manifests are now free
				Smatrix@labels[forWhich, forWhich] <- Smatrix@labels[fakeLatent,fakeLatent] # Give these manifests the fakeLatent's label				
			} else if (length(forWhich) > 1) {
				diag(Smatrix@values[forWhich, forWhich]) <- rep(Smatrix@values[fakeLatent,fakeLatent], length(forWhich)) # Move starting values
				diag(Smatrix@free[forWhich, forWhich]) <- rep(TRUE, length(forWhich)) # The manifests are now free
				# Give these manifests the fakeLatent's label
				diag(Smatrix@labels[forWhich, forWhich]) <- rep(Smatrix@labels[fakeLatent,fakeLatent], length(forWhich))
			}
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
		# fakeLatents <- array(0,0)
		latentVars <- realLatents
	}
	
	# -------------------------------------------------------------------------
	# IN DEVELOPMENT
	# Nonhomogenous error:
	# If Cerr is not proportional to I, a transformation is required before applying PPML
	# Cerr is built from the manifest variables in Smatrix and fakeLatents
	
	# Rinv will be applied later regardless -- By default, it is an identity matrix
	Rinv <- diag(1, length(manifestVars))
	
	# Assemble diagonal variances to check if NHEV transform is necessary
	diagVars <- rep(0,length(manifestVars))
	for (mvi in 1:length(manifestVars)) {
		if (Smatrix@free[mvi, mvi]) {
			diagVars[mvi] <- Smatrix@values[mvi, mvi]
		}
	}
	
	# for (fakeLatent in fakeLatents) {
		# diagVars[which(Amatrix@values[ ,fakeLatent] != 0)] <- Smatrix@values[fakeLatent, fakeLatent]
	# }
	
	# CHECKING: Any free values off the Cerr diagonal, or the Cerr diagonal is not homogeneous
	if ( any(as.logical(Smatrix@free[manifestVars,manifestVars] & !diag(rep(TRUE, length(manifestVars)))) ) ||
		(length(unique(diagVars)) > 1) ) {
		#(length(unique(diag(Smatrix@values[manifestVars, manifestVars]))) > 1) ) {
		
		bCERetVal <- buildCErr(Smatrix, manifestVars, constraints)
		if (is.null(bCERetVal))
			return(model)
		Cerr <- bCERetVal[[1]]
		relevantConstraints <- bCERetVal[[2]]
		
		# R <- t(chol(Cerr))  # OLD
		# R <- chol(Cerr)
		# Find transformation matrix R^-1
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
	
	# -------------------------------------------------------------------------
	# IN DEVELOPMENT: Data Missingness
	# If any data is missing, call function to split model in to submodels for
	# each pattern of missingness, transform each submodel, and create a
	# combined model.
	if (model@data@type == "raw" && any(is.na(model@data@observed))) {
		# In case an NHEV transformation has occurred, plug all the matrices
		# back in to a model to pass to the data missingness transform function
		MDmodel <- mxModel(model)
		MDmodel[[Aname]] <- Amatrix
		MDmodel[[Sname]] <- Smatrix
		MDmodel[[Fname]] <- Fmatrix
		if (!single.na(Mname)) {
			MDmodel[[Mname]] <- Mmatrix
		}
		MDmodel@constraints <- constraints
		if (!is.numeric(latentVars)) {
			MDmodel@latentVars <- latentVars
		} else if (length(latentVars) < (max(dim(Fmatrix@values)) - sum(Fmatrix@values)) ) {
			# if number of latentVars has been reduced, remove appropriate variable
			# names from the objective
			oldLatentVars <- NULL
			for (i in 1:max(dim(Fmatrix@values))) {
				if (!any(as.logical(Fmatrix@values[,i])))
					oldLatentVars <- c(latentVars, as.numeric(i))
			}
			MDmodel@objective@dims <- MDmodel@objective@dims[-setdiff(oldLatentVars, latentVars)]
		}
		
		
		MDmodel <- imxPPMLMissingData(MDmodel)
		if (!is.null(MDmodel))
		{
			return(MDmodel)
		} else {
			return(model)
		}	
	}

	###################
	#TRANSFORM SECTION#
	###################
	#get loadings from latents to manifest from A,S,F
	#transform A to E
	
	# Save original model for analytical solution
	# TODO: If filter deems that solution is only optimizable, don't do this
	originalModel <- model
	
	E <- solve(diag(nrow(Amatrix)) - Amatrix@values)
	#select only these columns of A which are real latents
	#get loadings matrix:= The part of A which goes from the latents to the manifests
	lambda <- as.matrix(E[manifestVars, realLatents])
	qrDecom <- qr(lambda)
	#check if loadings matrix is of full rank
	k <- length(realLatents)
	
	#if (qrDecom$rank != dim(lambda)[2] || k != qrDecom$rank){
	#	return(model)
	#}
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
	
	
	
	#transform data or cov
	if(model$data@type == "raw") {
		model$data@observed <- as.matrix(model@data@observed) %*% (Rinv) %*% t(Q)
		if (!is.numeric(manifestVars))
			colnames(model$data@observed) <- manifestVars
		else {
			colnames(model$data@observed) <- model@objective@dims[manifestVars]
		}
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
	# browser()
	# ANALYTICAL SOLUTION OR SPLIT MODEL
	
	# -------------------------------------------------------------------------
	# IN DEVELOPMENT -- ANALYTICAL SOLUTIONS
	# TODO: Test this with fakeLatents!
	# -------------------------------------------------------------------------
	# Only valid for saturated latent covariance matrix
	if ( solveType == "Solve" || solveType == "PartialSolve" ) {
		#|| ( all(diag(Smatrix@free[latentVars, latentVars])) && all( !(Smatrix@free[latentVars, latentVars] & diag(FALSE, nrow=length(latentVars), ncol=length(latentVars))) )) ) {
		
		# Error variance
		# For raw data
		if (model$data@type == 'raw') {
			# Calculate varE
			varE <- sum(model@data@observed[ ,(length(latentVars)+1):dim(model@data@observed)[[2]]]^2 )
			varE <- varE / ( (dim(model@data@observed)[[2]] - length(latentVars)) * dim(model@data@observed)[[1]] )
		}
		else if (model$data@type == 'cov') {
			varE <- sum(diag(model@data@observed[(length(latentVars)+1):dim(model@data@observed)[[1]], (length(latentVars)+1):dim(model@data@observed)[[2]]]))/(dim(model@data@observed)[[1]]-length(latentVars))
		}
		
		
		# Calculate lambdaInv
		lambdaInv <- solve(lambda[1:length(latentVars), ])
		
		# Analytical solution for latent means
		if (!is.null(Mmatrix)) {
			muLatent <- NULL
			if (model$data@type == 'raw') {
				# Calculate means of the latentVars
				# Covariance matrix of the first K variables in D'
				muLatent <- apply(as.matrix(model@data@observed[, 1:length(latentVars)]), 2, mean)
			} else if (model$data@type == 'cov') {
				muLatent <- model$data@means[1:length(latentVars)]
			}
			if ( is.numeric(latentVars) ) {
				Mmatrix@values[latentVars] <- lambdaInv %*% muLatent
			} else {
				# For some reason, indexing by character dimnames doesn't work
				# lapply in the index converts dimnames for latentvars in to appropriate numerical indices
				Mmatrix@values[ unlist(lapply(latentVars, function(x) { which(dimnames(Mmatrix)[[2]] == x) })) ] <- lambdaInv %*% muLatent
			}
			#Mmatrix@free[unlist(lapply(latentVars, function(x) { which(dimnames(Mmatrix)[[2]] == x) }))] <- FALSE
			
			if ( solveType == "Solve") {
				# Insert Mmatrix back in to original model
				originalModel[[Mname]] <- Mmatrix
			} else { # else if (solveType == "PartialSolve")
				# Analytical solution for the means:
				# Fix means values
				Mmatrix@free[unlist(lapply(latentVars, function(x) { which(dimnames(Mmatrix)[[2]] == x) }))] <- FALSE
				# Reinsert in to transformed model
				model[[Mname]] <- Mmatrix 
			}
		}
		if ( solveType == "Solve" ) {
			
			# Calculate covariance of the latentVars
			SigmaLatent <- NULL
			
			if (model$data@type == 'raw')
				SigmaLatent <- cov(as.matrix(model$data@observed[, 1:length(latentVars)]))
			else if (model$data@type == 'cov')
				SigmaLatent <- as.matrix(model$data@observed[1:length(latentVars), 1:length(latentVars)])
				
			# NOTE: Kludgy fix for symmetrization issues
			CLatent <- lambdaInv %*% (SigmaLatent - diag(varE, nrow=length(latentVars), ncol=length(latentVars))) %*% t(lambdaInv)
			Smatrix@values[latentVars, latentVars] <- (CLatent + t(CLatent))/2
			#Smatrix@free[latentVars, latentVars] <- FALSE
			
			# Insert in to S matrix
			Smatrix[manifestVars, manifestVars]@values <- diag(rep(varE, length(manifestVars))) # Implicitly: off-diagonal values of this submatrix are zeroed
			#diag(Smatrix[manifestVars, manifestVars]@free) <- rep(FALSE, length(manifestVars)) # Values are no longer free
			
			# browser()
			#  TODO: To make matrix spec work, shuffle the matrices so fake latents are the last variables in sequence,
			#  before variable indices are calculated
			#  Shuffle back at end, during restore
			# If there were fakelatents, fold fakelatents back out in the S matrix
			if (length(fakeLatents) > 0) {
				unfoldedS <- originalModel[[Sname]]
				unfoldedS@values[manifestVars, manifestVars] <- Smatrix@values[manifestVars, manifestVars]
				unfoldedS@values[latentVars, latentVars] <- Smatrix@values[latentVars, latentVars]
				for (fakeLatent in fakeLatents) {
					loads <- which(originalModel[[Aname]]@values[, fakeLatent])
					unfoldedS@values[fakeLatent, fakeLatent] <- Smatrix@values[loads, loads][1]
				}
			}
			
			# Reinsert Smatrix to model
			originalModel[[Sname]] <- Smatrix
			
			originalModel <- mxOption(originalModel, "UsePPML", "Solved")
			return(originalModel)
		} else { # (solveType == "PartialSolve")
			# Reinsert calculated varE in to Smatrix, fix values
			diag(Smatrix@values[manifestVars, manifestVars]) <- rep(varE, length(manifestVars))
			diag(Smatrix@free[manifestVars,manifestVars]) <- rep(FALSE, length(manifestVars))
			model[[Sname]] <- Smatrix
		}
	}
	
	# browser()
	
	###############
	#SPLIT SECTION#
	###############
	
	result <- NULL # Scope for result
	
	#Flatten spaces, commas out of model name to prevent problems
	model@name <- gsub(" ", "_", model@name)
	model@name <- gsub(",", "", model@name)
	
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
	leftmodel <- mxRename(model, paste(model@name,'_leftmodel',sep=""))	
	leftmodel <- selectSubModelFData(leftmodel, selectLatents, selectManifests)
	leftmodel <- mxOption(leftmodel, "UsePPML", "Split")
	
	# browser()
	
	if ( solveType == "PartialSolve" ) {
		result <- leftmodel
		if (model$data@type == 'raw') {
			result$data <- model$data
			result$data@observed <- result$data@observed[selectManifests]
			result$data@means <- result$data@means[selectManifests]
		}
		result <- mxOption(result, "UsePPML", "PartialSolved")
		
	} else if ( solveType == "Split") {
		
		# TODO: Rightmodel can be represented algebraically
		#rightmodel
		selectLatents <- latentVars[is.na(pmatch(latentVars, selectLatents))]
		selectManifests <- manifestVars[is.na(pmatch(manifestVars, selectManifests))]
		rightmodel <- mxRename(model, paste(model@name,'_rightmodel',sep=""))	
		rightmodel <- selectSubModelFData(rightmodel, selectLatents, selectManifests)
		rightmodel <- mxOption(rightmodel, "UsePPML", "Split")
	
		#reunite the two submodel
		result <-  mxModel(paste('PPMLTransformed_', model@name,sep=""), leftmodel, rightmodel)
	
	
		if (model$data@type == 'raw') {
			result$data <- model$data
		}
		
		modelnames <- c(paste(model@name,'_leftmodel',sep=""), paste(model@name,'_rightmodel',sep=""))
		objectives <- paste(modelnames, "objective", sep = ".")
		objectives <- paste(objectives, collapse = " + ")
		expression <- paste("mxAlgebra(", objectives, ", name = 'TotObj')", sep = "")
		algebra <- eval(parse(text=expression))
		objective <- mxAlgebraObjective("TotObj")
		result <- mxModel(result, objective, algebra)
		result <- mxOption(result, "UsePPML", "Split")
	} 
	
	return(result)
}

imxPPMLMissingData <- function(model) {
	# Need to name anonymous params to constrain across the submodels
	model <- omxNameAnonymousParameters(model)[[1]]
	
	Aname <- model@objective@A
	Sname <- model@objective@S
	Fname <- model@objective@F
	Mname <- model@objective@M
	
	# TODO / NOTE: Do we need to worry about screwy matrix-defined models? e.g., a case where
	# the latents and manifests are mixed in together in the matrices, rather than
	# the manifests and then the latents, in order
	
	# New patterns code
	# check through model@data@observed
	# Find unique missingness patterns
	patterns <- as.list(data.frame(t(unique(is.na(model@data@observed)))))
	patternLocs <- list()
	for ( patt in patterns) {
		patternLocs <- append(patternLocs, list(which(apply(is.na(model@data@observed), 1, function(row) { if (all(row == patt)) TRUE else FALSE } ))))
	}
	
#	# Check through model@data@observed, find the unique missingness patterns
#	# Automatically do first pattern
#	patterns <- list(is.na(model@data@observed[1, ]))
#	patternLocs <- list(1)
#	for ( r in 2:(dim(model@data@observed)[1]) ) {
#		thisPattern <- is.na(model@data@observed[r, ])
#		newPattern <- TRUE
#		for (p in 1:length(patterns)) {
#			if (all(patterns[[p]] == thisPattern)) {
#				patternLocs[[p]] <- c(patternLocs[[p]], r)
#				newPattern <- FALSE
#				break;
#			}
#		}
#		if (newPattern) {
#			patterns <- append(patterns, list(thisPattern))
#			patternLocs <- append(patternLocs, c(r))
#		}
#	}
	
	# Create a submodel for each missingness pattern with the appropriate
	# manifest variables removed
	bigModel <- mxModel(paste("(PPML Missing Data)", model@name))
	submodelNames <- c()
	PPMLAppliedCount <- 0
	for (m in 1:length(patterns)) {	
		newSubmodel <- mxModel(model)
		
		
		# Path specified versus matrix specified
		# Get label sets from data dimnames using LIVs
		removedMans <- (dimnames(newSubmodel@data@observed)[[2]])[patterns[[m]]]
		remainingMans <- (dimnames(newSubmodel@data@observed)[[2]])[!patterns[[m]]]
		remainingVars <- NULL # get remainingVars in scope
		remainingMansData <- remainingMans
		# NOTE: actually only need removedMans?
		if ( length(newSubmodel@manifestVars) == 0 && length(newSubmodel@latentVars) == 0 ) {
			remainingVars <- c(remainingMans, model@objective@dims[apply(newSubmodel$F@values, 2, sum) == 0])
			# Matrix specified - Use dims from objective to get location vector
			removedMansLoc <- c()
			for (removedMan in removedMans) {
				removedMansLoc <- c(removedMansLoc, which(model@objective@dims == removedMan))
			}
			removedMans <- sort(removedMansLoc)
			
			remainingMansLoc <- c()
			for (remainingMan in remainingMans) {
				remainingMansLoc <- c(remainingMansLoc, which(model@objective@dims == remainingMan))
			}
			remainingMans <- sort(remainingMansLoc)
			
			remainingVarsLoc <- c()
			for (remainingVar in remainingVars) {
				remainingVarsLoc <- c(remainingVarsLoc, which(model@objective@dims == remainingVar))
			}
			remainingVars <- sort(remainingVarsLoc)
		} else {
			# MAYBE:
			# NEED TO PUT REMAININGMANS AND REMAININGVARS IN THE PROPER ORDER
			# Is there ever a case where a path specified model produces matrices
			# with the manifests and latents mixed?
			remainingVars <- c(remainingMans, newSubmodel@latentVars)
		}
	
		# Fix matrices
		newSubmodel[[Aname]] <- newSubmodel[[Aname]][remainingVars, remainingVars]
		newSubmodel[[Sname]] <- newSubmodel[[Sname]][remainingVars, remainingVars]
		newSubmodel[[Fname]] <- newSubmodel[[Fname]][unlist(apply(newSubmodel[[Fname]]@values[ , remainingVars], 2, function(r) which(as.logical(r)))), remainingVars]
		if (!single.na(Mname)) {
			newSubmodel[[Mname]] <- newSubmodel[[Mname]][ , remainingVars]
		}
		
		# Remove in data, including means vectors
		# TODO / NOTE: Kludgy fix for the data matrix turning in to a vector when there's only one
		# manifest remaining -- should this actually happen?
		if (length(remainingMans) == 1 ) {
			newSubmodel@data@observed <- as.matrix(newSubmodel@data@observed[patternLocs[[m]], remainingMansData])
			colnames(newSubmodel@data@observed) <- remainingMansData
		} else {
			# NORMAL CASE -- should probably be this all the time
			newSubmodel@data@observed <- newSubmodel@data@observed[patternLocs[[m]], remainingMansData]
		}
		
		if (!single.na(newSubmodel@data@means)) {
			newSubmodel@data@means <- newSubmodel@data@means[, remainingMans]
		}
		
		# Fix list of manifestVars
		if ( !(length(newSubmodel@manifestVars) == 0 && length(newSubmodel@latentVars) == 0) ) {
			newSubmodel@manifestVars <- remainingMans
		} else {
			newSubmodel@objective@dims <- newSubmodel@objective@dims[-removedMans]
		}
		
		# Rename each submodel
		newSubmodel@name <- paste('PPML_MG_Submodel', m, sep="")
		
		# Transform each of these submodels with PPML
		newSubmodel <- imxTransformModelPPML(newSubmodel)
		submodelNames <- c(submodelNames, newSubmodel@name)
		
		
		# Check if PPML was applicable on the submodel
		if (!is.null(newSubmodel@options$UsePPML) && newSubmodel@options$UsePPML == "Applied") { # TODO : Fix UsePPML check here
			PPMLAppliedCount <- PPMLAppliedCount + 1
		}
		
		# Combine these submodels in to a larger model
		bigModel <- mxModel(bigModel, newSubmodel)
	}
	
	# Need to check if any of the submodels were successfully transformed; if not, reject transformation
	if (PPMLAppliedCount == 0) {
		return(NULL)
	}
	
	# Objective = sum
	objectives <- paste(submodelNames, "objective", sep = ".")
	objectives <- paste(objectives, collapse = " + ")
	expression <- paste("mxAlgebra(", objectives, ", name = 'SumObjective')", sep = "")
	algebra <- eval(parse(text=expression))
	objective <- mxAlgebraObjective("SumObjective")
	bigModel <- mxModel(bigModel, objective, algebra)
	bigModel <- mxOption(bigModel, "UsePPML", "Applied")
	
	# Return this model
	return(bigModel)
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
	if (!(is.numeric(selectManifests) && is.numeric(selectLatents))) {
		submodel$F <- model$F[selectManifests, Aindices]
	} else {
		# Matrix-specified:
		# for each Aindice along the columns, find the row that contains a 1
		# keep the columns corresponding to Aindices and those rows
		Frows <- c()
		for (i in 1:dim(model$F)[[2]]) {
			if (any(Aindices == i))
				Frows <- c(Frows, which(as.logical(model$F@values[ ,i])))
		}
		submodel$F <- model$F[Frows, Aindices]
	}
	
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
	
		# Restore dimnames to the new covariance matrix
		if (!(is.numeric(selectManifests) && is.numeric(selectLatents))) {
			# path specified
			colnames(submodel$data@observed) <- selectManifests
			rownames(submodel$data@observed) <- selectManifests
		} else {
			# matrix specified
			colnames(submodel$data@observed) <- colnames(model$data@observed)[selectManifests]
			rownames(submodel$data@observed) <- rownames(model$data@observed)[selectManifests]
		}	
		
		# If means data exists, pull out the appropriate manifest variables
		if (!single.na(submodel@data@means)){
			submodel@data@means <- t(as.matrix(submodel@data@means[1,selectManifests]))
		}
	}
	return(submodel)
}

buildCErr <- function(Smatrix, manifestVars, constraints) {
	# BASIC: If there are no constraints, then the model is not valid for PPML
	if (length(constraints) == 0)
		return(NULL)

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
		# NOTE: Constraint checking algorithm should only return(NULL) if the
		# constraint involves at least one relevant label. Otherwise, it can 
		# just ignore the constraint.
		for (term in 1:length(constraint@formula[[2]])) {
			# If the term is a label, do checks
			if (is.character(constraint@formula[[2]][[term]]) ) {
				if (any(as.character(constraint@formula[[2]][[term]]) == errLabels)) {
					if (hasRelevantLabelLeft || hasIrrelevantLabelLeft) {
						# Only one label per side allowed
						return(NULL)
					}
					hasRelevantLabelLeft <- TRUE
				} else {
					hasIrrelevantLabelLeft <- TRUE
					if (hasRelevantLabelLeft) {
						# Can't have irrelevant and relevant labels in the same constraint
						return(NULL)
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
						return(NULL)
					}
					if (hasIrrelevantLabelRight || hasRelevantLabelRight) {
						# Can't have more than one label on a side
						return(NULL)
					}
					hasRelevantLabelRight <- TRUE
				} else {
					hasIrrelevantLabelRight <- TRUE
					if (hasRelevantLabelLeft) {
						# Can't have a relevant label left and an irrelevant right
						return(NULL)
					}
					if (hasRelevantLabelRight) {
						# Can't have two labels on one side
						return(NULL)
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
			return(NULL)
		}
		
		# Left and right hand sides must be of length 3 or 1
		# NOTE: Not a very valuable test, necessary at all?
		# NOTE: Could fold in to below loop
		if (!( (length(constraint@formula[[2]]) == 1 || length(constraint@formula[[2]]) == 3) &&
			(length(constraint@formula[[3]]) == 1 || length(constraint@formula[[3]]) == 3) ) ) {
			return(NULL)
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
					return(NULL)
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
					return(NULL)
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
								return(NULL)
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
								return(NULL)
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
								return(NULL)
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
								return(NULL)
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
		return(NULL)
	}
	
	# Create Cerr
	for (errLabel in errLabels) {
		Cerr[Smatrix@labels[manifestVars, manifestVars] == errLabel] <- errValues[which(errLabels == errLabel)]
	}
	
	return(list(Cerr, relevantConstraints))
}

imxRestoreResultPPML <- function(model, result) {
	omxCheckTrue( (length(grep("PPML", result@name, fixed = TRUE)) > 0) )
	
	manifestVars <- NULL # Scope
	latentVars <- NULL
	if (length(model@manifestVars > 0)) {	
		manifestVars <- model@manifestVars
		latentVars <- model@latentVars
	} else {
		manifestVars <- c()
		for (i in 1:max(dim(model$F@values))) {
			if (any(as.logical(model$F@values[,i])))
				manifestVars <- c(manifestVars, as.numeric(i))
			else
				latentVars <- c(latentVars, as.numeric(i))
		}
	}
	
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
		for (latentVar in latentVars) {
			if (sum(model@matrices$A@values[,latentVar]) == 1) {
				if (sum(model@matrices$A@values[manifestVars, latentVar]) == 1) {
					fakeLatents <- c(fakeLatents, latentVar)
				}
			}
		}
		
		# Deal with fakelatents:
		if (length(fakeLatents) > 0) {
			# Rebuild the model, folding fakeLatents in to the S matrix
			for (fakeLatent in fakeLatents) {
				forWhich <- which(model$A@values[ ,fakeLatent] != 0)
				Smatrix@values[forWhich, forWhich] <- Smatrix@values[fakeLatent,fakeLatent]
				Smatrix@free[forWhich, forWhich] <- TRUE
				Smatrix@labels[forWhich, forWhich] <- Smatrix@labels[fakeLatent,fakeLatent]
			}
			# Remove fakeLatents from A, S, and F matrices
			remainingVars <- c(manifestVars, setdiff(latentVars, fakeLatents))
			Smatrix <- Smatrix[remainingVars, remainingVars]
		}
		

		Cerr <- buildCErr(Smatrix, manifestVars, model@constraints)
		
		errScale <- paramValues[which(paramLabels == "_PPML_NHEV_ErrParam")]
		
		newParamLabels <- array(0,0)
		newParamValues <- array(0,0)
		for (errParam in setdiff(unique(model@matrices$S@labels), c(paramLabels, NA))) {
			newParamLabels <- c(newParamLabels, errParam)
			newParamValues <- c(newParamValues, errScale * setdiff(unique(Cerr[[1]][which(Smatrix@labels[manifestVars, manifestVars] == errParam)]), NA))
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

imxPPML <- function(model, flag) {
	if (!(isS4(model) && is(model, "MxModel"))) {
		stop("Argument 'model' must be MxModel object")
	}
	if (!(length(flag) == 1 && is.logical(flag) && !is.na(flag))) {
		stop("Argument 'flag' must be TRUE or FALSE")
	}
	# Objective exists / must be a RAM objective
	objective <- model$objective
	if(is.null(objective) || !is(objective, "MxRAMObjective")) {
		return(NA)
	}
	if (flag) {
		model@options$UsePPML <- "Yes"
	} else {
		model@options$UsePPML <- "No"
	}
	return(model)
}

single.na <- function(a) {
    return((length(a) == 1) &&
        (is.list(a) || is.vector(a) || is.matrix(a)) &&
        (is.na(a) == TRUE))
}	

# DEVELOPMENT FUNCTIONS

# Function could be extended to be more general, but for now is just for comparing fits for models,
# their PPML transforms, and the corresponding reversed PPML transforms
# checkLL : Disable for NHEV models, log likelihood will not match
# checkByName : Disable for some matrix models (some models specified with latents and manifests mixed together)
#    NOTE: Must name all parameters in these models for the check to be complete
imxCheckFitsPPML <- function(res1, res2, checkHessians = TRUE, checkLL, checkByName) {
	# Check -2logLLs versus each other
	if (checkLL)
		omxCheckCloseEnough(res1@output$Minus2LogLikelihood, res2@output$Minus2LogLikelihood, 0.01)
		
	if (checkByName) {
		
		for (label in labels(res1@output$estimate)) {
			omxCheckCloseEnough(res1@output$estimate[label], res2@output$estimate[label], 0.01)
		}
		
		oneInTwo <- lapply(labels(res1@output$estimate), function(l) { which(labels(res2@output$estimate) == l) })
		
		omxCheckCloseEnough(res1@output$standardErrors, as.matrix(res2@output$standardErrors[unlist(oneInTwo)]), .01)
		
	} else {
		# Check parameters versus each other
		# Parameters will be in the same order, so a simple iteration should work
		for ( i in 1:length(res1@output$estimate) ) {
			omxCheckCloseEnough(res1@output$estimate[i], res2@output$estimate[i], 0.01)
		}
		
		# Check standard errors versus each other
		# Should just be able to check the two vectors directly against each other
		omxCheckCloseEnough(res1@output$standardErrors, res2@output$standardErrors, 0.01)
	}

	
	# TODO: Check if @expMeans are replicated...
	# 

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

imxTestPPML <- function(model, checkLL = TRUE, checkByName = FALSE) {
	res1 <- mxRun(imxPPML(model, FALSE), suppressWarnings = TRUE) # Standard fit
	#res2 <- mxRun(imxTransformModelPPML(model))	# PPML fit
	res2 <- mxRun(imxPPML(model, TRUE), suppressWarnings = TRUE)	# PPML fit
	# browser()
	#omxCheckTrue( (length(grep("PPML", res2@name, fixed = TRUE)) > 0) )
	#omxCheckTrue( is(res2@objective, "MxAlgebraObjective") ) # redundant, could comment this out
	# First model not transformed
	omxCheckTrue( !(!is.null(res1@options$UsePPML) && (res1@options$UsePPML == "Solved" || res1@options$UsePPML == "PartialSolved")) )
	# Second model transformed
	omxCheckTrue( (!is.null(res2@options$UsePPML) && (res2@options$UsePPML == "Solved" || res2@options$UsePPML == "PartialSolved")) )
	#res3 <- imxRestoreResultPPML(model, res2)	# Reverse transform PPML
	
	# NOTE: Not checking Hessians anymore, they're never very close -- not sure that
	# a comparison is actually meaningful
	# NOTE: checkLL parameter can turn off log likelihood checking for this test,
	# useful for non-homogeneous error variance cases where the transformed log 
	# likelihood is different.
	imxCheckFitsPPML(res1, res2, checkLL = checkLL, checkByName = checkByName, checkHessians = FALSE) # Check standard fit vs PPML model
	
	# NOTE Always check likelihood for this test, to make sure everything was plugged in correctly
	#imxCheckFitsPPML(res1, res3, checkLL = TRUE, checkByName = checkByName, checkHessians = FALSE)	# Check standard fit vs reverse transformed PPML model
	
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
