	# Cases: Solvable, and Optimizable

	# Solvable Case
	# I -- Name anonymous parameters in model
	# II -- Apply pretransforms to model
	# III -- Apply transform to model
	# IV --- Solve using transformed model
	# V --- Plug parameters back in to named model
	# VI --- Plug values matrix from named model back in to original model
	# return original model

	# Optimizable case - Partial Solution
	# I -- Name anonymous parameters in model
	# II -- Apply pretransforms to model
	# III --- Apply transform to model
	# IV -- Split model, solve for variance, plug back in to leftmodel
	# V -- Call optimizer on model
	# VI -- Extract parameters and reinsert to named model
	# VII -- Plug values matrix from named model back in to original model
	# return original model

	# Optimizable case - Missing Data
	# I -- Name anonymous parameters in model
	# II -- Apply pretransforms to model
	# III -- Split model for missingness patterns
	# IV -- Apply transform to each applicable submodel
	# V -- Pass overall model to optimizer
	# VI -- Extract parameters and reinsert to named model
	# VII -- Plug values matrix from named model back in to original model

	# Optimizable case - Split Model
	# I -- Name anonymous parameters in model
	# II -- Apply pretransforms to model
	# III -- Apply transform to model
	# IV -- Split model
	# V -- Call optimizer on model
	# VI -- Extract parameters and reinsert to named model
	# VII - Plug values matrix from named model back in to original model

	# SO:
	# I -- Name anonymous parameters in model
	# II -- Apply pretransforms to model
	# III --  Missing data split, transform each applicable submodel -or-
	# 				Apply transform to model
	# IV -- Solve using transformed model -or-
	#       solve for variance using transformed model, optimize -or-
	#				split model -or-
	#				nothing if missing data split already happened
	#	V -- If partial solve or split (MD or otherwise) call optimizer
	# VI -- Extract parameters and reinsert to named model
	# VII - Plug values matrix from named model back in to original model


setClass("PPMLSolveType", 
	representation(
		result = "character",
		isMultiLayer = "logical",
		isMatrixSpecified = "logical",
		hasNHEV = "logical",
		Cerr = "matrix",
		relevantConstraints = "vector",
		hasMissingness = "logical",
		hasNoLatentCovariances = "logical",
		latentCovNotSaturated = "logical",
		hasFixedExpectedLatentMean = "logical", 
		manifestVars = "vector",
		latentVars = "vector",
		fakeLatents = "vector"
	),
	prototype(
		result="Check",
		isMultiLayer = FALSE,
		isMatrixSpecified = FALSE,
		hasNHEV = FALSE,
		Cerr = NULL,
		relevantConstraints = NULL,
		hasMissingness=FALSE,
		hasNoLatentCovariances = FALSE,
		latentCovNotSaturated = FALSE,
		hasFixedExpectedLatentMean = FALSE,
		manifestVars = character(0),
		latentVars = character(0),
		fakeLatents = character(0)
	)
)

##' imxPPML
##'
##' Potentially enable the PPML optimization for the given model.
##' 
##' @param model the MxModel to evaluate
##' @param flag whether to potentially enable PPML
imxPPML <- function(model, flag=TRUE) {
	if (!(isS4(model) && is(model, "MxModel"))) {
		stop("Argument 'model' must be MxModel object")
	}
	if (!(length(flag) == 1 && is.logical(flag) && !is.na(flag))) {
		stop("Argument 'flag' must be TRUE or FALSE")
	}
	# Expectation exists / must be a RAM expectation
	expectation <- model$expectation
	if(is.null(expectation) || !is(expectation, "MxExpectationRAM")) {
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

PPMLTransformModel <- function(model.original) {
	# Name anonymous parameters in model
	pair <- (omxNameAnonymousParameters(model.original))
	model.named <- pair[[1]]

	solveType <- PPML.CheckApplicable(model.named)


	# PPML not applicable to model
	if (single.na(solveType))
	{
		model.original@options$UsePPML <- "Inapplicable"
		return(model.original)
	}

	model.pretransformed <- model.named

	# Apply pretransforms to model

	# Remove changes made to data by the mxRun function
	if (!is.null(model.pretransformed$data@means))
		model.pretransformed@data <- mxData(numObs=model.pretransformed@data@numObs, observed=model.pretransformed@data@observed, type=model.pretransformed@data@type, means=as.vector(model.pretransformed@data@means))
		

	
	# Fold in fakelatents
	if (length(solveType@fakeLatents))
	{
		model.pretransformed <- PPML.Pre.FakeLatents(model.pretransformed, solveType)
		solveType@latentVars <- setdiff(solveType@latentVars, solveType@fakeLatents)
		# If matrix specified, specified by indices -- need to refigure manifestVars and latentVars in solveType
		if ( is.numeric(solveType@latentVars) && is.numeric(solveType@manifestVars) )
		{
			# Iterate through F matrix to find which columns are manifests (1 present) and
			# which columns are latents (no 1 present in column)
			Fmatrix <- model.pretransformed[[model.pretransformed$expectation@F]]
			manifestVars <- numeric(0)
			latentVars <- numeric(0)
			for (i in 1:dim(Fmatrix@values)[[2]]) {
				if (any(as.logical(Fmatrix@values[,i])))
					manifestVars <- c(manifestVars, as.numeric(i))
				else
					latentVars <- c(latentVars, as.numeric(i))	
			}
			solveType@manifestVars <- manifestVars
			solveType@latentVars <- latentVars
		}
	}	


	model.noFakeLatents <- model.pretransformed
	# Build Cerr (with nonfakelatented model)
	# Also detects if the model is not a valid PPML NHEV model
	# Aborts PPML transform if so

	solveType <- PPML.Pre.UpdateSolveType(model.noFakeLatents, solveType)

	if (single.na(solveType))
	{
		model.original@options$UsePPML <- "Inapplicable"
		return(model.original)
	}
	

	# Fix NHEV
	if (solveType@hasNHEV)
	{

		# Fix up the model
		model.pretransformed <- PPML.Pre.FixNHEV(model.pretransformed, solveType)
	}


	# Apply transforms to model
	if (solveType@hasMissingness)
	{
		# Missing Data Case
		model.transformed <- PPMLMissingData(model.pretransformed, solveType)

		# Check to make sure PPML is applicable to model (shouldn't ever happen here)
		if (single.na(model.transformed))
		{
			model.original@options$UsePPML <- "Inapplicable"
			return(model.original)
		}

		model.transformed@options$UsePPML <- "No"
		results <- mxRun(model.transformed)
	}
	else
	{
		pair <- PPML.Transform(model.pretransformed, solveType)

		model.transformed <- pair[[1]]
		lambda <- pair[[2]]

		if (solveType@result == "PartialSolve" || solveType@result == "Solve")
			model.transformed <- PPML.SolveOrPartialSolve(model.transformed, lambda, solveType)

		if (solveType@result == "PartialSolve" )
		{

			model.transformed@options$UsePPML <- "No"
			# Optimize partially solved
			results <- mxRun(model.transformed)

			# Restore free values
			results <- omxSetParameters(results, labels=names(omxGetParameters(model.named)), free=TRUE)
		}
		else if (solveType@result == "Split")
		{
			# Split the model if necessary
			model.transformed <- PPML.Split(model.transformed, solveType)
			model.transformed@options$UsePPML <- "No"

			# Optimize the split model
			results <- mxRun(model.transformed)
		}
		else
			results <- model.transformed
	}

	if (solveType@hasNHEV)
	{
		results <- PPML.Post.UnfoldNHEV(results, model.noFakeLatents, solveType)
	}

	



	# Extract parameters and reinsert to named model
	params <- omxGetParameters(results)
	model.named <- omxSetParameters(model.named, names(params), values=as.vector(params))
	
	
	# Plug values matrix from named model back in to original model
	model.original[[model.original$expectation@S]]@values <- model.named[[model.original$expectation@S]]@values
	if (!single.na(model.original$expectation@M))
		model.original[[model.original$expectation@M]]@values <- model.named[[model.original$expectation@M]]@values

	# Preserve output data and indicate what type of solution was used
#	model.original@output <- results@output
	model.original@options$UsePPML <- solveType@result

	# Return original model with solution
	return(model.original)
}




PPML.CheckApplicable <- function(model) {
	solveType = new("PPMLSolveType")	


	# Explicitly, don't use PPML -OR- no explicit instruction to use PPML and default is no
	if (model@options$UsePPML == "No" || (is.null(model@options$UsePPML) && getOption("mxOptions")$UsePPML == "No"))
	{
		#print("PPML abort: disabled")
		return(NA)
	}
 
	#is the model a RAM model?
	# NOTE: This could be made redundant by implementing the transform call
	# from the MxExpectationRAM function
	expectation <- model$expectation
	if(is.null(expectation) || !is(expectation, "MxExpectationRAM")) {
		#print("PPML abort: Not a RAM model")
		return(NA)
	}

	# Extract RAM matrices
	Aname <- expectation@A
	Sname <- expectation@S
	Fname <- expectation@F
	Mname <- expectation@M
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
		#print("PPML abort: Matrices not MxMatrices")
		return(NA)
	}
	
	# Make sure data is present, and raw or covariance
	if( !is.null(model$data) ) {
		if ( !(model$data@type == "raw" || model$data@type == "cov") )  {
			#print("PPML abort: Need raw or covariance matrix data")
			return(NA)
		}
	} else {
		#print("PPML abort: Model has no data")
		return(NA)
	}
	
	# Structure matrix must be fixed
	# TODO: Optimization possible for cases where only part of the structure matrix is fixed
	if(any(Amatrix@free)) {
		#print("PPML abort: Structure matrix not fixed")
		return(NA)	
	}


	# If the model is matrix specified, uses numeric indices instead of dimnames
	# Then for split models: splits the expectation function in the split section
	manifestVars <- NULL
	latentVars <- NULL

	solveType@isMatrixSpecified <- length(model@manifestVars) == 0 && length(model@latentVars) == 0
	if (solveType@isMatrixSpecified) {

		if (length(model$expectation@dims) == 1 && single.na(unique(model$expectation@dims)))	# no dims anywhere NOTE: is this necessary?
		{
			#print("PPML abort: Model is missing dims")
			return(NA)
		}

		# Iterate through F matrix to find which columns are manifests (1 present) and
		# which columns are latents (no 1 present in column)
		for (i in 1:dim(Fmatrix@values)[[2]]) {
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

	# Put manifestVars, latentVars in to solveType
	solveType@manifestVars <- manifestVars
	solveType@latentVars <- latentVars

	### Call latent classifier function to classify latents
	classifiedLatents <- PPMLClassifyLatents(Amatrix, Smatrix, latentVars)
	if (single.na(classifiedLatents))
	{
		#print("PPML abort: Latents of improper form")
		return(NA)
	}
	realLatents <- classifiedLatents[[1]]
	fakeLatents <- classifiedLatents[[2]]
	solveType@fakeLatents <- fakeLatents
	rootLatents <- classifiedLatents[[3]]
	nonRootLatents <- classifiedLatents[[4]]


	# Check if expected means vector for latents are FREE
	# If all are not free, then solution doesn't apply
	# Use solution as best guess, but don't fix them for partial solution
	# If all are not free, analytical solution becomes partial solution
	if (!single.na(Mname))
	{
		if (!all(Mmatrix@free[,realLatents]))
			solveType@hasFixedExpectedLatentMean <- TRUE
	}
	
	# DEVELOPMENT -- Some of these models should be transformable
	# TODO: Multilayered models
	if (length(nonRootLatents) > 0) {
		solveType@isMultiLayer <- TRUE
	}
	
	# (I-Lambda)^-1 -- Structure matrix must be invertible
	# Check for multi-layered models
	# PROFILING: Expensive
#	if ( det( diag(dim(Amatrix@values)[[1]]) - Amatrix@values ) == 0 ) {
#		return(NA)
#	}
	
	# BASIC: Model must have manifestvars and latentvars and more manifests than latents
	if (length(realLatents) == 0 || length(manifestVars) == 0 || 
		length(realLatents) >= length(manifestVars)) {
		#print("PPML abort: Must have manifests, latents, and more manifests than latents")
		return(NA)
	}
	
	# DATA MISSINGNESS CHECK	
	if (model$data@type == "raw" && any(is.na(model$data@observed))) {
		# Check missingness patterns: if none of the patterns have more manifests than latents,
		# then PPML cannot be applied at all
		if ( (length(manifestVars) - min(apply(is.na(model$data@observed),1,sum))) <= length(realLatents) ) {
			#print("PPML abort: No missingness patterns are transformable")
			return(NA)
		}
		
		# Might be solvable, but for now, split only
		solveType@hasMissingness <- TRUE
		#solveType <- "Split" 
	}

	#  Partial Solve possible when all covariances between latents are fixed to zero but all
	# variances of the latents are free
	#  Analytical solution still possible for cases with one latent
	#  All covariances are fixed, zero
	# Conditional below:
	#  More than one latent
	#  There are no free covariances between the latents
	hasNoFreeCovariances <- !any(Smatrix@free[realLatents, realLatents] & !diag(TRUE, length(realLatents), length(realLatents)) )
	#  All latent covariances are (effectively) zero
	allLatentCovsZero <- all( abs(Smatrix@values[realLatents, realLatents] - diag(diag(Smatrix@values[realLatents,realLatents]))) < .001 ) 

	if ( hasNoFreeCovariances && allLatentCovsZero ) {

		if ( all(as.logical(diag(Smatrix@free[realLatents, realLatents]))) ) {
			# All latent variances are free 
			#if (solveType != "Split") solveType <- "PartialSolve"
			solveType@hasNoLatentCovariances <- TRUE
		} else {
				# One or more latent variance is fixed
				#print("PPML abort: One or more latent variance is fixed")
				#return(NA)
				solveType@latentCovNotSaturated <- TRUE
				
		}

	} else if ( !all(Smatrix@free[realLatents, realLatents]) ) {
#		# If some less structured combination of variances and covariances of the latents
#		# is fixed, probably not applicable 
#		print("PPML abort: Inapplicable variance/covariance structure (Fixedness)")
#		return(NA)
		solveType@latentCovNotSaturated <- TRUE
	}

	# Multilayer doesn't work with missingness
	if (solveType@isMultiLayer && solveType@hasMissingness)
		return(NA)

	# Determine solveType results
	if (!solveType@hasMissingness) {
		# No data missingness
		if ((solveType@hasNoLatentCovariances || solveType@latentCovNotSaturated) && (length(latentVars) - length(fakeLatents)) > 1) {
			solveType@result <- "Split"
		} else {
			if (solveType@hasFixedExpectedLatentMean)
				solveType@result <- "PartialSolve"
			else
				solveType@result <- "Solve"
		}
	} else {
		# Data missingness
		solveType@result <- "MissingData"
	}
	
	# Although not explicit in the code,
	# solveType <- "Solve" happens implicitly when Smatrix@free[latentVars, latentVars] is all free
	# (saturated covariance matrix)
	
	# BROKEN CASE: Covariance data + means, not analytical solution
	# BLOCK THIS CASE
	# Covariance data
	# Means data present
	# Solvetype is not "Solve" == analytical soln
	if (model$data@type == "cov" && !single.na(model$data@means) && solveType@result != "Solve")
		return(NA)

	return(solveType)
	
}

PPML.Check.UseOptimizer <- function(opt)
{
		if (is.null(opt))
			return(TRUE)
		if (!(opt == "No" || opt == "Inapplicable"))
			return(FALSE)
		return(TRUE)
}


PPML.Pre.UpdateSolveType <- function(model, solveType)
{
	manifestVars <- solveType@manifestVars
	latentVars <- solveType@latentVars
	
	# Extract RAM matrices
	expectation <- model$expectation
	Aname <- expectation@A
	Sname <- expectation@S
	Fname <- expectation@F
	Mname <- expectation@M
	Amatrix <- model[[Aname]]
	Smatrix <- model[[Sname]]
	Fmatrix <- model[[Fname]]
	Mmatrix <- NULL
	if (!single.na(Mname)) {
		Mmatrix <- model[[Mname]]
	}


	maybeNHEV <- FALSE

	# NHEV CHECK -- PART I
	# Error matrix := Smatrix[manifestVars, manifestVars]
	# With any variances in fake latents in the appropriate diagonal spots

	# TODO: In NHEV cases, constraints that aren't part of the error covariance matrix -- can PPML
	# be applied with any other constraint on the model, or will further constraints always break it?
	
	# Any free values off the diagonal in the manifest covariance matrix?
	hasFreeValuesOffDiagonal <- any(as.logical( Smatrix@free[manifestVars,manifestVars] & !diag(rep(TRUE, length(manifestVars))) ) )

	# If there are unlabeled params in the error matrix, not even an NHEV case
	# -- NHEV algorithm relies on constraints between labeled params
	# Can't do anything with this
	hasUnlabeledFreeParams <- !all(Smatrix@free[manifestVars,manifestVars] == !is.na(Smatrix@labels[manifestVars,manifestVars]))
	if (hasUnlabeledFreeParams) # IMPLICIT: >1 manifestVar, otherwise PPML not applicable anyways. SO the two must be constrained together somehow
	{
#		print("PPML abort: Unlabeled free error parameters")
		return(NA)
	}

	# Need constraints for NHEV, but can also break model if not relevant to error matrix
	hasConstraints <- length(model@constraints) > 0

	# Figure out if there's more than one label on the diagonals of the error matrix

	hasMultipleDiagErrorLabels <- length(unique(diag(Smatrix@labels[manifestVars, manifestVars]))) > 1

	# Actual checking
	# First case, non-diagonal matrix
	if (hasFreeValuesOffDiagonal) {
		# Need constraints for the error matrix to vary as matrixC*scalarVarE
		if (!hasConstraints)
		{
#			print("PPML abort: Unconstrained error variance")
			return(NA)
		}
		
		if (!hasMultipleDiagErrorLabels)
		{
			# Matrix diagonals are all equal
			# Need for the free values off the diagonal to be labeled differently
			# Otherwise you get Identity + additional 1s scattered in the matrix,
			# which are not (in general?) invertible
			diagLabels <- unique(diag(Smatrix@labels[manifestVars,manifestVars]))
			allLabels <- unique(as.vector(Smatrix@labels[manifestVars,manifestVars]))

			hasDifferentLabelsOffDiag <- as.logical(length(setdiff(allLabels, c(diagLabels,NA))))
			if (!hasDifferentLabelsOffDiag)
			{
#				print("PPML abort: Constraints make non-positive-definite error matrix")
				return(NA)
			}
		}

		maybeNHEV <- TRUE # Actually, right now, just "Maybe"
		# Second half of check happens after fakeLatents are folded in to model,
		# Partway through the pretransform phase
		# TODO: Maybe make a second parameter, maybeNHEV, for clarity?
	} 
	# Second case, diagonal matrix
	else
	{
		if (hasMultipleDiagErrorLabels)
		{
			# Need constraints for the error matrix to vary as matrixC*scalarVarE
			if (!hasConstraints)
			{
#				print("PPML abort: Constraints required on error matrix for NHEV transform")
				return(NA)			
			}
			maybeNHEV <- TRUE
		}
		# IMPLICIT else: Diagonal matrix where all error variances are constrained
		# to equality --> standard case, not NHEV

	}


	if (!maybeNHEV)
		return(solveType)


	hasIrrelevantConstraints <- FALSE

	# NHEV Check -- Part II
	# Generates the Cerr matrix
	# Leave this for latest possible in case the model is found to be PPML-inapplicable
	# This is possibly really slow, and should be avoided if at all possible

	constraints <- model@constraints

	# IMPORTANT: NHEV and Missingness algorithms do not seem to be compatible
	if (solveType@hasMissingness)
		return(NA)

	# Call buildCErr to build CErr matrix
	bCERetVal <- buildCErr(model[[model$expectation@S]], solveType@manifestVars, constraints)
	if (is.null(bCERetVal))
		return(NA)
	Cerr <- bCERetVal[[1]]
	relevantConstraints <- bCERetVal[[2]]

	hasIrrelevantConstraints <- length(relevantConstraints) < length(constraints)
	# DEVELOPMENT
	# If there are irrelevant constraints, analytical solution is probably wrong
	# Might as well try it?
	# TODO: See if this breaks any models
	# /DEVELOPMENT
	if (hasIrrelevantConstraints)
		solveType@result <- "Split" # Downgrade solveType to split
	# IMPLICIT ELSE: leaves solveType as it was before (NHEV will never improve it)

	# buildCErr returns null if there's an invalid constraint
	if (is.null(Cerr))
		return(NA)

	# Noninvertible Cerr
	# --> Not a valid PPML model, need to be able to invert Cerr
	if (det(Cerr) == 0)
		return(NA)

	solveType@Cerr <- Cerr
	solveType@relevantConstraints <- relevantConstraints
	solveType@hasNHEV <- TRUE

	return(solveType)

}








PPMLClassifyLatents <- function(Amatrix, Smatrix, latentVars) {
	# Classify latents
	fakeLatents <- array(0,0)
	rootLatents <- array(0,0)
	for(latentVar in latentVars){
		# Get nonzero loadings from the asymmetric matrix
		loads <- which(Amatrix@values[,latentVar] != 0)
		
		# Fake latents only load on to one manifest variable
		if ( length(loads) == 1) {
			# Load must have value 1
			# Otherwise, folding fakelatent back out becomes difficult and less graceful methods are required
			if (Amatrix@values[loads[1], latentVar] == 1) {
				# All loadings from the latent must be on manifestVars
				#  NOTE: Latents that covary with other latents but have no loadings on to anything?
				#  Lambda is probably singular in these cases.
				if (length(which(Amatrix@values[latentVars, latentVar] != 0)) == 0) { # No loadings on latentVars
					# Latent must not covary with any manifests or other latents
					if ( Smatrix@free[latentVar,latentVar] && length(which(Smatrix@free[,latentVar])) == 1) {
						# All of these manifestVars must have no inherent error variances /
						# loads are all on to manifests without their own error variances
						if ( !Smatrix@free[loads, loads] && Smatrix@values[loads, loads] == 0 ) {
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

	return(list(realLatents, fakeLatents, rootLatents, nonRootLatents))
}




PPML.Pre.FakeLatents <- function(model, solveType)
{
	Aname <- model$expectation@A
	Sname <- model$expectation@S
	Fname <- model$expectation@F
	Mname <- model$expectation@M

	Amatrix <- model[[Aname]]
	Smatrix <- model[[Sname]]
	Fmatrix <- model[[Fname]]
	Mmatrix <- NULL
	if (!single.na(Mname)) {
		Mmatrix <- model[[Mname]]
	}

	latentVars <- solveType@latentVars	
	manifestVars <- solveType@manifestVars
	fakeLatents <- solveType@fakeLatents

	# Rebuild the model, folding fakeLatents in to the S matrix
	for (fakeLatent in fakeLatents) {
		forWhich <- which(Amatrix@values[ ,fakeLatent] != 0) # Manifest to which the fakeLatent corresponds
		loadVal <- Amatrix@values[forWhich]

		Smatrix@values[forWhich, forWhich] <- Smatrix@values[fakeLatent,fakeLatent] # Move starting value
		Smatrix@free[forWhich, forWhich] <- TRUE # The manifest is now free
		Smatrix@labels[forWhich, forWhich] <- Smatrix@labels[fakeLatent,fakeLatent] # Give the manifest the fakeLatent's label				

	}

	# Remove fakeLatents from A, S, and F matrices
	remainingVars <- c(manifestVars, setdiff(latentVars, fakeLatents))
	
	# Get remainingVars in correct order
	if (solveType@isMatrixSpecified)
	{
		remainingVars <- sort(remainingVars)
	}
	else
	{
		sortedVars <- c()
		vars <- colnames(Fmatrix)
		for (var in vars)
		{
			if (any(remainingVars == var))
				sortedVars <- c(sortedVars, var)
		}
		remainingVars <- sortedVars
	}
	
	Amatrix <- Amatrix[remainingVars, remainingVars]
	Smatrix <- Smatrix[remainingVars, remainingVars]
	Fmatrix <- Fmatrix[ ,remainingVars]
	if (!single.na(Mname)) {
		Mmatrix <- Mmatrix[,remainingVars]
	}

	model[[Aname]] <- Amatrix
	model[[Sname]] <- Smatrix
	model[[Fname]] <- Fmatrix
	if (!single.na(Mname))
		model[[Mname]] <- Mmatrix

	# Fix expectation function/latent variable list
	if (!is.numeric(latentVars)) { # Path-specified
		model@latentVars <- setdiff(latentVars, fakeLatents)
	} else { # Matrix-specified
		# Repair dims
		model$expectation@dims <- model$expectation@dims[-fakeLatents]
	}

	return(model)

}




PPML.Pre.FixNHEV <- function(model, solveType) {
	Aname <- model$expectation@A
	Sname <- model$expectation@S
	Mname <- model$expectation@M

	Amatrix <- model[[Aname]]
	Smatrix <- model[[Sname]]
	Mmatrix <- NULL
	if (!single.na(Mname)) {
		Mmatrix <- model[[Mname]]
	}

	constraints <- model@constraints
	manifestVars <- solveType@manifestVars
	latentVars <- solveType@latentVars

	Cerr <- solveType@Cerr
	relevantConstraints <- solveType@relevantConstraints


	# Find transformation matrix R^-1
	Rinv <- (solve(chol(Cerr)))
	
	# Apply it to original Cerror to get transformed C'error
	# (Actually, just generate appropriately sized identity matrix)
	CerrPrime <- diag(rep(1, dim(Cerr)[1]))
	
	Smatrix@values[manifestVars, manifestVars] <- CerrPrime # Reinsert to Smatrix
	
	# Set all labels for the errors to the same label
	Smatrix@free[manifestVars, manifestVars] <- FALSE
	Smatrix@labels[manifestVars, manifestVars] <- NA
	for (manifestVar in manifestVars) {
		if (Smatrix@values[manifestVar, manifestVar] != 0) {
			Smatrix@labels[manifestVar, manifestVar] <- '_PPML_NHEV_ErrParam' #unique(clabels)
			Smatrix@free[manifestVar, manifestVar] <- TRUE # NOTE: ? maybe
		}
	}

	# Transform structure matrix
	Amatrix[manifestVars, latentVars]@values <- t(Rinv) %*% Amatrix[manifestVars, latentVars]@values

	# Reinsert transformed matrices to model
	model[[Aname]] <- Amatrix
	model[[Sname]] <- Smatrix

	# TODO: Is there a missing transform on the M matrix?
	if (!single.na(Mname))
		model[[Mname]] <- Mmatrix

	# Remove relevant constraints
	for (relevantConstraint in relevantConstraints) {
		constraints <- setdiff(constraints, list(model@constraints[[relevantConstraint]]))
	}

	model@constraints <- constraints

	# Transform data
	if (model$data@type == "raw")
	{
		# Transform raw data
		model$data@observed <- as.matrix(model$data@observed) %*% (Rinv)
		if (!is.numeric(manifestVars))
			colnames(model$data@observed) <- manifestVars
		else {
			colnames(model$data@observed) <- model$expectation@dims[manifestVars]
		}

	}
	else if (model$data@type == "cov")
	{
		# Transform variance data
		model$data@observed <- t(Rinv) %*% as.matrix(model$data@observed) %*% (Rinv)

		if (!single.na(model$data@means))
			model$data@means <- as.matrix(model$data@means) %*% (Rinv)
		
		# Restore dimnames
		if (!(is.numeric(manifestVars) && is.numeric(latentVars))) {
			# For path-specified models:
			colnames(model$data@observed) <- manifestVars
			rownames(model$data@observed) <- manifestVars
			if (!single.na(model$data@means))
				colnames(model$data@means)    <- manifestVars
		} else {
			# For matrix-specified models
			colnames(model$data@observed) <- model$expectation@dims[manifestVars]
			rownames(model$data@observed) <- model$expectation@dims[manifestVars]
			if (!single.na(model$data@means))
				colnames(model$data@means)    <- model$expectation@dims[manifestVars]
		}
	}
	return(model)
}




buildCErr <- function(Smatrix, manifestVars, constraints) {
	# BASIC: If there are no constraints, then the model is not valid for PPML
	if (length(constraints) == 0)
		return(NULL)

	# Build Cerr

	# Get labels of manifestVars
	errLabels <- setdiff(unique(as.vector(Smatrix@labels[manifestVars,manifestVars])), NA)
	# Generate appropriately sized and dimnamed matrix
	Cerr <- matrix(data=0, nrow=length(manifestVars), ncol=length(manifestVars), dimnames=list(manifestVars, manifestVars))
	
	# Label[1] * Factor == Label[2]
	conLabels <- as.character(array(0,0))
	conFactors <- as.numeric(array(0,0))
	
	relevantConstraints <- as.character(array(0,0))
	
	for (constraint in constraints) {
		# NOTE: Necessary to get a length check on the constraint?
		
		# Constraint is relevant if left and right hand side each involve
		# one label from the error matrix, and no other labels
		# If a constraint is a relationship between terms in the error matrix
		# and terms outside of the error matrix, then the model is not of the
		# correct form.
		# NOTE: This would not be true in the case of a label on a fixed value
		# TODO: Maybe propagate fixedness out before applying NHEV transform?
		
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

PPML.Post.UnfoldNHEV <- function(result, model.noFakeLatents, solveType)
{
	Sname <- model.noFakeLatents$expectation@S
	manifestVars <- solveType@manifestVars

	# Extract solved error param
	errParam <- omxGetParameters(result)[['_PPML_NHEV_ErrParam']]
	
	# Scale Cerr matrix appropriately
	scaledCerr <- solveType@Cerr * errParam

	# Pull out labels, free associated with Cerr matrix
	labelMatrix <- model.noFakeLatents[[Sname]]@labels[manifestVars,manifestVars]
	freeMatrix <- model.noFakeLatents[[Sname]]@free[manifestVars,manifestVars]

	# Put labels, free, scaled Cerr back in to result
	result[[Sname]]@values[manifestVars,manifestVars] <- scaledCerr
	result[[Sname]]@labels[manifestVars,manifestVars] <- labelMatrix
	result[[Sname]]@free[manifestVars,manifestVars] <- freeMatrix

	return(result)
	
}


PPMLMissingData <- function(model, solveType) {
	# Need to name anonymous params to constrain across the submodels
	Aname <- model$expectation@A
	Sname <- model$expectation@S
	Fname <- model$expectation@F
	Mname <- model$expectation@M



	latMeanVarNames <- character(0)
	if (!single.na(Mname))
	{
		latMeanVarNames <- model[[Mname]]@labels[,solveType@latentVars]
		# Strip out NA (occurs when some latent means are fixed, not free anyways)
		latMeanVarNames <- latMeanVarNames[ which(!is.na(latMeanVarNames)) ]
	}

	errVarName <- diag(model[[Sname]]@labels[solveType@manifestVars, solveType@manifestVars])[[1]]

	# Check through model$data@observed
	# Find unique missingness patterns
	patterns <- as.list(data.frame(t(unique(is.na(model$data@observed)))))
	patternLocs <- list()
	patternCounts <- list()
	for ( patt in patterns) {
		patternLocs <- append(patternLocs, list(which(apply(is.na(model$data@observed), 1, function(row) { if (all(row == patt)) TRUE else FALSE } ))))
	}
	
	# Create a submodel for each missingness pattern with the appropriate
	# manifest variables removed
	bigModel <- mxModel(paste("(PPML Missing Data)", model@name))
	submodelNames <- c()
	PPMLAppliedCount <- 0
	PPMLSolvedCount <- 0

	meanVarE <- 0
	meanLatentMeans <- numeric(0)
	if (!single.na(Mname))
		meanLatentMeans <- vector(mode="numeric", length=length(latMeanVarNames))

	for (m in 1:length(patterns)) {	
		newSubmodel <- mxModel(model)
		newSolveType <- solveType
		
		# Path specified versus matrix specified
		# Get label sets from data dimnames using LIVs
		removedMans <- (dimnames(newSubmodel$data@observed)[[2]])[patterns[[m]]]

		remainingMans <- (dimnames(newSubmodel$data@observed)[[2]])[!patterns[[m]]]
		remainingVars <- NULL # get remainingVars in scope
		remainingMansData <- remainingMans
		# NOTE: actually only need removedMans?
		if ( solveType@isMatrixSpecified ) {
			remainingVars <- c(remainingMans, model$expectation@dims[solveType@latentVars])
			# Matrix specified - Use dims from expectation to get location vector
			removedMansLoc <- c()
			for (removedMan in removedMans) {
				removedMansLoc <- c(removedMansLoc, which(model$expectation@dims == removedMan))
			}
			if (length(removedMans)) removedMans <- sort(removedMansLoc)
			
			remainingMansLoc <- c()
			for (remainingMan in remainingMans) {
				remainingMansLoc <- c(remainingMansLoc, which(model$expectation@dims == remainingMan))
			}
			remainingMans <- sort(remainingMansLoc)
			
			remainingVarsLoc <- c()
			for (remainingVar in remainingVars) {
				remainingVarsLoc <- c(remainingVarsLoc, which(model$expectation@dims == remainingVar))
			}
			remainingVars <- sort(remainingVarsLoc)
		} else {

			remainingVars <- c(remainingMans, newSolveType@latentVars)
		}
	
		# If any manifests are missing from the data, remove those manifests from the submodel
		if (!is.null(removedMans))
		{
			# Fix matrices
			newSubmodel[[Aname]] <- newSubmodel[[Aname]][remainingVars, remainingVars]
			newSubmodel[[Sname]] <- newSubmodel[[Sname]][remainingVars, remainingVars]
			newSubmodel[[Fname]] <- newSubmodel[[Fname]][unlist(apply(newSubmodel[[Fname]]@values[ , remainingVars], 2, function(r) which(as.logical(r)))), remainingVars]
			if (!single.na(Mname)) {
				newSubmodel[[Mname]] <- newSubmodel[[Mname]][ , remainingVars]
			}

			if (!single.na(newSubmodel$data@means)) {
				newSubmodel$data@means <- newSubmodel$data@means[, remainingMans]
			}
		
			# Fix list of manifestVars
			if ( !solveType@isMatrixSpecified ) {
				# PATH
				newSubmodel@manifestVars <- remainingMansData
				newSolveType@manifestVars <- remainingMansData
			} else {
				# MATRIX
				dimIndices <- sort(c(remainingMans, solveType@latentVars))
				newSubmodel$expectation@dims <- newSubmodel$expectation@dims[dimIndices]
				newSolveType@manifestVars <- remainingMans
				# Need to fix up indices, as items have potentially been removed, shifting indices
				for (remVar in rev(sort(removedMans)))
				{
					# Decrement indices > remVar by one
					newSolveType@manifestVars <- unlist(lapply(newSolveType@manifestVars, function (x) { if (x > remVar) return(x-1) else return(x) } ))
					newSolveType@latentVars <- unlist(lapply(newSolveType@latentVars, function (x) { if (x > remVar) return(x-1) else return(x) } ))
				}

			}
		}

		# Remove pattern-appropriate manifests in data, including means vectors
		# NOTE: Kludgy fix for the data matrix turning in to a vector when there's only one
		# manifest remaining -- should this actually happen?
		if (length(remainingMans) == 1 ) {
			newSubmodel$data@observed <- as.matrix(newSubmodel$data@observed[patternLocs[[m]], remainingMansData])
			newSubmodel$data@numObs <- dim(newSubmodel$data@observed)[[1]]
			colnames(newSubmodel$data@observed) <- remainingMansData
		} else {
			# NORMAL CASE -- should probably be this all the time
			newSubmodel$data@observed <- newSubmodel$data@observed[patternLocs[[m]], remainingMansData]
			newSubmodel$data@numObs <- dim(newSubmodel$data@observed)[[1]]
			colnames(newSubmodel$data@observed) <- remainingMansData
		}

		
		# Rename each submodel
		newSubmodel@name <- paste('PPML_MG_Submodel', m, sep="")
		
		# Transform each of these submodels with PPML,
		# if applicable
		if (length(newSolveType@manifestVars) > length(newSolveType@latentVars))
		{
			pair <- PPML.Transform(newSubmodel, newSolveType)
			newSubmodel <- PPML.SolveOrPartialSolve(pair[[1]], pair[[2]], newSolveType)
			newSubmodel <- PPML.Split(newSubmodel, newSolveType)
		}		
		submodelNames <- c(submodelNames, newSubmodel@name)
		
		# Check if PPML was applicable on the submodel
		if (!is.null(newSubmodel@options$UsePPML) && 
			(newSubmodel@options$UsePPML == "Solved" || newSubmodel@options$UsePPML == "PartialSolved" || newSubmodel@options$UsePPML == "Split")) { 
			PPMLAppliedCount <- PPMLAppliedCount + 1

			# Below should be zero, as solution should not be applied to missing data submodels
			if (newSubmodel@options$UsePPML == "Solved") 
				PPMLSolvedCount <- PPMLSolvedCount + 1

			# Track mean of the solved error variances to produce best estimate for starting point
			meanVarE <- meanVarE + omxGetParameters(newSubmodel)[errVarName]
			# Track mean of the solved latent means to produce best estimate for starting point
			if (!single.na(Mname))
			{
				latMeans <- omxGetParameters(newSubmodel)[latMeanVarNames]
				for (i in 1:length(latMeanVarNames))
				{
					latName <- latMeanVarNames[i]
					if (!single.na(meanLatentMeans[latName]))
						meanLatentMeans[i] <- meanLatentMeans[i] + meanLatentMeans[latName]
				}
				
				meanLatentMeans <- meanLatentMeans
			}

		} 

		# Set UsePPML to "No" to allow optimization to occur
		newSubmodel@options$UsePPML <- "No"
		
		# Combine these submodels in to a larger model
		bigModel <- mxModel(bigModel, newSubmodel)
	}
	
	# Need to check if any of the submodels were successfully transformed; if not, reject transformation
	if (PPMLAppliedCount == 0) {
		return(NA)
	}

	meanVarE <- meanVarE/PPMLAppliedCount
	meanLatentMeans <- meanLatentMeans/PPMLAppliedCount
	bigModel <- omxSetParameters(bigModel, c(errVarName, latMeanVarNames), values=c(meanVarE,meanLatentMeans))
	
	# Objective = sum
	# TODO: Combine algebra objectives from PPML-transformed models?
	#objectives <- paste(submodelNames, "objective", sep = ".")
	objectives <- paste(submodelNames, "fitfunction", sep = ".")
	objectives <- paste(objectives, collapse = " + ")
	expression <- paste("mxAlgebra(", objectives, ", name = 'SumObjective')", sep = "")
	algebra <- eval(parse(text=expression))
	objective <- mxFitFunctionAlgebra("SumObjective")
	bigModel <- mxModel(bigModel, objective, algebra)
	bigModel <- mxOption(bigModel, "UsePPML", "MissingData")


	return(bigModel)
}





PPML.Transform <- function(model, solveType)
{
	Aname <- model$expectation@A
	Sname <- model$expectation@S
	Fname <- model$expectation@F
	Mname <- model$expectation@M

	Amatrix <- model[[Aname]]
	Smatrix <- model[[Sname]]
	Fmatrix <- model[[Fname]]
	Mmatrix <- NULL
	if (!single.na(Mname)) {
		Mmatrix <- model[[Mname]]
	}

	manifestVars <- solveType@manifestVars
	latentVars <- solveType@latentVars

	#transform A to E	
	E <- solve(diag(nrow(Amatrix)) - Amatrix@values)
	#select only these columns of A which are real latents
	#get loadings matrix:= The part of A which goes from the latents to the manifests
	lambda <- as.matrix(E[manifestVars, latentVars])
	
	qrDecom <- qr(lambda)
	
	#check if loadings matrix is of full rank
	k <- length(latentVars)
	#if (qrDecom$rank != dim(lambda)[2] || k != qrDecom$rank){
	#	return(model)
	#}
	#last check passed, can start modifying model

	#calculate upper triangle matrix
	#orthogonal
	Q <- t(qr.Q(qrDecom, complete = TRUE)) 

	#transform loadings matrix
	lambda <- Q %*% lambda
	#set all rows from k+1 to zero
	lambda[(k + 1) : nrow(lambda), ] <- 0
	#insert new loadings matrix in A
	Amatrix@values[manifestVars,latentVars] <- lambda
	Amatrix@values[latentVars,latentVars] <- 0 # For latents predicting latents


	# Reinsert transformed matrices back in to model
	model[[Aname]] <- Amatrix
	model[[Sname]] <- Smatrix
	model[[Fname]] <- Fmatrix
	if (!single.na(Mname)) {
		model[[Mname]] <- Mmatrix
	}
	
	#transform data or cov
	if(model$data@type == "raw") {
		model$data@observed <- as.matrix(model$data@observed) %*% t(Q)

		# Restore dimnames
		if (!is.numeric(manifestVars))
			colnames(model$data@observed) <- manifestVars
		else {
			colnames(model$data@observed) <- model$expectation@dims[manifestVars]
		}
	}
	else if(model$data@type == "cov") {
		# Transform variance data
		model$data@observed <- (Q) %*% as.matrix(model$data@observed) %*% t(Q)

		# Restore dimnames
		if (!solveType@isMatrixSpecified) {
			# For path-specified models:
			colnames(model$data@observed) <- manifestVars
			rownames(model$data@observed) <- manifestVars
		} else {
			# For matrix-specified models
			colnames(model$data@observed) <- model@expectation@dims[manifestVars]
			rownames(model$data@observed) <- model@expectation@dims[manifestVars]
		}

		# Transform means data, if it exists
		if (!single.na(model$data@means)){
			model$data@means <- as.matrix(model$data@means) %*% t(Q)
			# Restore dimnames
			if (!solveType@isMatrixSpecified)
				colnames(model$data@means) <- manifestVars	# Path-specified
			else
				colnames(model$data@means) <- model$expectation@dims[manifestVars] # Matrix-specified
		}
		
	}

	return(list(model, lambda))

}


PPML.SolveOrPartialSolve <- function(model, lambda, solveType)
{
	# Extract matrices from model
	expectation <- model$expectation

	Aname <- expectation@A
	Sname <- expectation@S
	Fname <- expectation@F
	Mname <- expectation@M
	Amatrix <- model[[Aname]]
	Smatrix <- model[[Sname]]
	Fmatrix <- model[[Fname]]
	Mmatrix <- NULL
	if (!single.na(Mname)) {
		Mmatrix <- model[[Mname]]
	}

	# Make option "Solve" or "PartialSolve"
	model@options$UsePPML <- solveType@result

	manifestVars <- solveType@manifestVars
	latentVars <- solveType@latentVars


	# Analytical solution for error variance
	# For raw data
	if (model$data@type == 'raw') {
		# Calculate varE
		K <- length(latentVars)
		M <- dim(model$data@observed)[[2]]
		N <- dim(model$data@observed)[[1]]
		varE <- sum(model$data@observed[ ,(K+1):M]^2 )
		varE <- varE / ( (M - K) * N )
	}
	else if (model$data@type == 'cov')
	{
		# Two components to the error variance in covariance data models:
		# The (post-transform) variances of the variables in the error model...
		K <- length(latentVars)
		M <- dim(model$data@observed)[[1]]
		varE <- sum(diag(as.matrix(model$data@observed[(K+1):M, (K+1):M]))) / (M-K)
		# And the sum of the square of the (post-transform) means of those variables
		if (!single.na(model$data@means))
		{
			# Need correction of N/(N-1) -- trying to generate sample variance, not population
			# For the variances pulled out of the transformed data cov matrix, this has been accounted for
			# However, denominator of means is N, not N-1
			N <- model$data@numObs
			varE <- varE + sum(model$data@means[,(K+1):M]^2) / (M-K) * N/(N-1)
		}
	}
	

	# Calculate lambdaInv
	lambdaInv <- solve(lambda[1:length(latentVars), ])
	
	# Analytical solution for means of latents
	if (!is.null(Mmatrix)) {
		muLatent <- NULL
		# Extract untransformed means
		if (model$data@type == 'raw') {
			# Calculate means of the manifest vars
			muLatent <- apply(as.matrix(model$data@observed[, 1:length(latentVars)]), 2, mean)
		} else if (model$data@type == 'cov') {
			# Manifest means are already specified
			muLatent <- model$data@means[1:length(latentVars)]
		}

		# Apply transform to means
		solvedM <- lambdaInv %*% muLatent
		if ( solveType@isMatrixSpecified ) {

			# If any of the latent means are fixed, can't just blindly put in the solved values
			# So, check:
			if (solveType@hasFixedExpectedLatentMean)
			{
				# Don't overwrite value of fixed latent
				for (i in 1:length(solveType@latentVars))
				{
					latentVar <- solveType@latentVars[i]
					if (Mmatrix@free[latentVar])	
						Mmatrix@values[latentVar] <- solvedM[i]
				}
			}
			else
				Mmatrix@values[latentVars] <- solvedM

		} else {
			# For some reason, indexing by character dimnames doesn't work here
			# fix: lapply call in the index converts dimnames for latentvars in to appropriate numerical indices
			latentIndices <- unlist(lapply(latentVars, function(x) { which(dimnames(Mmatrix)[[2]] == x) }))

			# If any of the latent means are fixed, can't just blindly put in the solved values
			# So, check:
			if (solveType@hasFixedExpectedLatentMean)
			{
				# If some of them are fixed, iterate through, checking each one
				# If it's not fixed, insert the solved value as a best guess
				for (i in 1:length(latentIndices))
				{
					latentIndex <- latentIndices[i]
					if (Mmatrix@free[latentIndex])
						Mmatrix@values[latentIndex] <- solvedM[i]
				}
			}
			else
				Mmatrix@values[ latentIndices ] <- lambdaInv %*% muLatent

		}
		

		# Analytical solution for the means:
		# Fix means values so optimizer doesn't bother with them (if partialsolved; if free, don't bother)
		# -- Don't fix these if, in the original model, any of the expected latent means are not free to vary
		# -- Don't fix in MissingData case, as each solved submodel will have a different solution
		# 	-- Just plug them in, use the average as a best guess, and optimize simultaneously
		if (solveType@result == "PartialSolve" && !solveType@hasFixedExpectedLatentMean)
			Mmatrix@free[unlist(lapply(latentVars, function(x) { which(dimnames(Mmatrix)[[2]] == x) }))] <- FALSE

		# Insert Mmatrix back in to transformed model
		model[[Mname]] <- Mmatrix
	}


	if (solveType@result == "PartialSolve" || solveType@result == "MissingData")
	{
		# Reinsert calculated varE in to Smatrix, fix values
		diag(Smatrix@values[manifestVars, manifestVars]) <- rep(varE, length(manifestVars))

		if (solveType@result == "PartialSolve")
		# If the model is being partially solved and then optimized, fix the error variances
			diag(Smatrix@free[manifestVars,manifestVars]) <- rep(FALSE, length(manifestVars))

		# Reinsert Smatrix in to model
		model[[Sname]] <- Smatrix

		# SO: Returns model with solution to error variances, latent means
		# This model must then be optimized
		return(model)
	}


	# IMPLICIT  solveType@result == "Solve"
	# Full analytical solution is possible
		
	# Calculate covariance of the latentVars
	SigmaLatent <- NULL
	
	if (model$data@type == 'raw')
	{
		# Using R's covariance matrix function does not work
		# SigmaLatent <- cov(as.matrix(model$data@observed[, 1:length(latentVars)]))

		# Must make covariance matrix of the manifests "by hand"
		RowxRowT <- apply(as.matrix(model$data@observed[,1:length(latentVars)]), 1, function(row) { as.matrix(row) %*% t(as.matrix(row)) })
		if (!is.matrix(RowxRowT))
			sigma <- sum(RowxRowT)
		else
			sigma <- matrix(apply(as.matrix(RowxRowT), 1, sum), nrow=length(latentVars), ncol=length(latentVars))
		sigma <- sigma/dim(model$data@observed)[[1]]

		# Covariance matrix in SigmaLatent
		SigmaLatent <- sigma - as.matrix(muLatent) %*% t(as.matrix(muLatent))

		
	}
	else if (model$data@type == 'cov')
		SigmaLatent <- as.matrix(model$data@observed[1:length(latentVars), 1:length(latentVars)])
		
	# Calculate latent covariance matrix and reinsert in to the S matrix
	# NOTE: Kludgy fix for symmetrization issues
	CLatent <- lambdaInv %*% (SigmaLatent - diag(x=varE, nrow=length(latentVars), ncol=length(latentVars))) %*% t(lambdaInv)
	Smatrix@values[latentVars, latentVars] <- (CLatent + t(CLatent))/2		
	
	# Reinsert error variance in to S matrix
	Smatrix@values[manifestVars, manifestVars] <- diag(x=varE, nrow=length(manifestVars), ncol=length(manifestVars)) # Implicitly: off-diagonal values of this submatrix are zeroed
	
	# Reinsert Smatrix to model
	model[[Sname]] <- Smatrix

	# Returns fully solved model
	return(model)
}


PPML.Split <- function(model, solveType) {
	# Extract matrices from model
	expectation <- model$expectation
	
	Aname <- expectation@A
	Sname <- expectation@S
	Fname <- expectation@F
	Mname <- expectation@M
	Amatrix <- model[[Aname]]
	Smatrix <- model[[Sname]]
	Fmatrix <- model[[Fname]]
	Mmatrix <- NULL
	if (!single.na(Mname)) {
		Mmatrix <- model[[Mname]]
	}

	manifestVars <- solveType@manifestVars
	latentVars <- solveType@latentVars

	#Flatten spaces, commas out of model name to prevent problems with expectation functions
	model@name <- gsub(" ", "_", model@name)
	model@name <- gsub(",", "", model@name)
	
	#leftmodel
	#calculate A realLatents indices for left model
	#first k manifest vars
	selectManifests <- manifestVars[1:length(latentVars)]
	if (solveType@isMatrixSpecified)
		selectData <- model$expectation@dims[selectManifests]
	else
		selectData <- selectManifests

	selectLatents <- latentVars


	leftmodel <- selectSubModelFData(model, selectLatents, selectManifests)
	leftmodel <- mxRename(leftmodel, paste(leftmodel@name, "_leftmodel", sep=""))
	leftmodel@options$UsePPML <- "No"

	# DEVELOPMENT: This will cause problems for models where the left and right model are solved simultaneously
	# Will need to put in a conditional	
	# In cases where the left and right model are solved simultaneously,
	# Need to wrap both in a bigModel and put the data in to that bigModel
	# TODO: Move this block to selectsubmodelfdata ?
	if (model$data@type == 'raw') { # RAW DATA
		# Extract model
		leftmodel$data <- model$data
		# Select relevant data
		leftmodel$data@observed <- as.matrix(leftmodel$data@observed[,selectData])
		if (!single.na(leftmodel$data@means))
			leftmodel$data@means <- leftmodel$data@means[selectData]

		# Restore dimnames
		colnames(leftmodel$data@observed) <- selectData
		if (!single.na(leftmodel$data@means))
			colnames(leftmodel$data@means) <- selectData
	}



	#rightmodel
	selectLatents <- latentVars[is.na(pmatch(latentVars, selectLatents))]
	selectManifests <- manifestVars[(length(latentVars)+1):length(manifestVars)]

	if (solveType@isMatrixSpecified)
		selectData <- model$expectation@dims[selectManifests]
	else
		selectData <- selectManifests
	
	if (model$data@type == 'raw') # RAW DATA
	{
		errData <- sum(model$data@observed[ ,selectData]^2)

		numErrData <- dim(model$data@observed)[[1]] * length(selectData)
	}
	else	# COV DATA
	{
		N <- model$data@numObs
		errData <- sum(diag(as.matrix(model$data@observed))[selectData])

		if (!single.na(model$data@means))
		{
			# Correction factor of N/N-1 at the end is because we're trying
			# to find a sample variance, not a population variance
			# The denominator of the mean is N, not N-1, so this must be corrected
			# to get the sample variance.
			# TRIED: No factor -- diff 9.194952e-2
			# N/(N-1) factor -- diff -8.420196e-2
			# N/(N-1) factor inside ()^2
			# (N-1)/N factor
			# sqrt(N/(N-1)) factor
			# (N+1)/N
			# (N-1)/(N-2)
					
			errData <- errData + sum(model$data@means[,selectData]^2)
		}

		errData <- errData * N
		numErrData <- N * length(selectData)
	}
	
	# IN DEVELOPMENT: Rightmodel represented algebraically
	# DEV: Old code, useful for left and right model solved simultaneously?
		# rightmodel <- mxRename(model, paste(model@name,'_rightmodel',sep=""))	
		# rightmodel <- selectSubModelFData(rightmodel, selectLatents, selectManifests)
		# rightmodel <- mxOption(rightmodel, "UsePPML", "Split")
	
	
	### Create separate 1x1 "error matrix" to be referenced by expectation, constrain to error variance with a label
	# Find error label for Smatrix
	errorLabel <- unique(diag(Smatrix@labels[manifestVars, manifestVars])) # NOTE: Works as long as there's only one error label
	# Create mxMatrix for error parameter
	#  (Value will be first appropriately labeled parameter in the Smatrix, although all should be set the same anyways)
	varE <- mxMatrix( "Full", nrow=1, ncol=1, labels=errorLabel, values=leftmodel$S@values[[ which(leftmodel$S@labels == errorLabel)[[1]] ]], free=TRUE, name="varE" )
	
	
	### Create algebra expectation
	sumSqData <- mxMatrix("Full", nrow=1, ncol=1, values = errData, name = "sumSqData")
	numData <- mxMatrix("Full", nrow=1, ncol=1, values = numErrData, name="numData")
	
	# (Build this algebra:)
	# 		errorLLAlgebra <- mxAlgebra(numData * log(eval(errorLabel)) + sumSqData/eval(errorLabel), name="errorLL")
	# Build algebra string
	expression <- paste("mxAlgebra(numData[1,1] * log(varE[1,1]) + sumSqData[1,1]/varE[1,1], name='errorLL')", sep="")
	# Create algebra
	errorLLAlgebra <- eval(parse(text=expression))
	

	# DEV: Old code, useful for simultaneously solving left and right models?
		#reunite the two submodel
		# result <-  mxModel(paste('PPMLTransformed_', model@name,sep=""), leftmodel, rightmodel)
	
	
	# Build algebra string
	#expression <- paste(paste(model@name,'_leftmodel',sep=""), "objective", sep=".") # "model@name_leftmodel.objective"
	expression <- paste(paste(model@name,'_leftmodel',sep=""), "fitfunction", sep=".") # "model@name_leftmodel.fitfunction"
	expression <- paste(c(expression, "errorLL"), collapse = " + ") # "model@name_leftmodel.objective + errorLL"
	expression <- paste("mxAlgebra(", expression, ", name='TotObj')", sep="") # mxAlgebra(model@name_leftmodel.objective + errorLL, name='TotObj')
	# Create algebra
	objAlgebra <- eval(parse(text=expression)) 


	#algObjective <- mxAlgebraObjective("TotObj") # Create algebraobjective that call this algebra
	#result <- mxModel(name=model@name, submodels=leftmodel, objective=algObjective, objAlgebra, sumSqData, numData, errorLLAlgebra, varE)
	algFit <- mxFitFunctionAlgebra("TotObj") # Create algebraobjective that call this algebra
	#result <- mxModel(name=model@name, submodels=leftmodel, objective=algObjective, objAlgebra, sumSqData, numData, errorLLAlgebra, varE)
	
	
	result <- mxModel(name=model@name, submodels=leftmodel, fitfunction=algFit, objAlgebra, sumSqData, numData, errorLLAlgebra, varE)
	

#		DEV: Will be necessary for models where left and right model are solved simultaneously,
#		and right model cannot be represented algebraically		
#		# Include data in combination model for raw data
#		if (model$data@type == 'raw') {
#			result$data <- model$data
#		}
	
	result <- mxOption(result, "UsePPML", "Split")
	
	# DEV: old code -- could be useful for true splitmodel case, where error and
	# leftmodel need to be optimized simultaneously
		# modelnames <- c(paste(model@name,'_leftmodel',sep=""), paste(model@name,'_rightmodel',sep=""))
		# objectives <- paste(modelnames, "objective", sep = ".")
		# objectives <- paste(objectives, collapse = " + ")
		# expression <- paste("mxAlgebra(", objectives, ", name = 'TotObj')", sep = "")
		# algebra <- eval(parse(text=expression))
		# objective <- mxAlgebraObjective("TotObj")
		# result <- mxModel(result, objective, algebra)
	
	return(result)
}



#@author:	Julian Karch
#		jk3nq@virginia.edu
#@idea>		Timo von Oertzen
#		timo@virginia.edu
selectSubModelFData <- function(model, selectLatents, selectManifests) {
	Aindices <- append(selectManifests, selectLatents)
	if (is.numeric(Aindices))
		Aindices <- sort(Aindices)

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
		submodel$expectation@dims <- submodel$expectation@dims[Aindices]
	}
	
	if(!is.null(submodel$M)){
		submodel$M <- model$M[1,Aindices]
		# submodel$M@values <- t(submodel$M@values)
		# submodel$M@labels <- t(submodel$M@labels)
		# submodel$M@free <- t(submodel$M@free)
		# submodel$M@lbound <- t(submodel$M@lbound)
		# submodel$M@ubound <- t(submodel$M@ubound)
	}
	
	if (is.numeric(selectManifests))
		selectData <- model$expectation@dims[selectManifests]
	else
		selectData <- selectManifests

	if(model$data@type == "raw"){
		submodel$data <- NULL
#		submodel$data@observed <- as.matrix(submodel$data@observed[,selectManifests])
	} else if (model$data@type == "cov") {
		# Pull out the proper manifest variables from the covariance matrix
		# to create the appropriate covariance matrix for the submodel
		submodel$data@observed <- as.matrix(submodel$data@observed[selectData,selectData])

		# If means data exists, pull out the appropriate manifest variables
		if (!single.na(submodel$data@means))
			submodel$data@means <- t(as.matrix(submodel$data@means[1,selectData]))
	
		# Restore dimnames to the new covariance matrix
		colnames(submodel$data@observed) <- selectData
		rownames(submodel$data@observed) <- selectData
		if (!single.na(submodel$data@means))
			colnames(submodel$data@means) <- selectData
			
	}
	return(submodel)
}


###############################################################################
### TESTING PPML
###############################################################################

PPML.Tool.CheckPPMLDidEqualOrBetter <- function(omxRes, ppmlRes, tolerance)
{
	#diff <- omxRes@output$Minus2LogLikelihood - ppmlRes@output$Minus2LogLikelihood
	diff <- omxRes@output$minimum - ppmlRes@output$minimum

	if (diff < -tolerance)
		#stop("ppmlRes@output$Minus2LogLikelihood not less than or equal to omxRes@output$Minus2LogLikelihood within tolerance of ", tolerance,"\n",omxRes@output$Minus2LogLikelihood," vs ",ppmlRes@output$Minus2LogLikelihood)	
		stop("ppmlRes@output$minimum not less than or equal to omxRes@output$minimum within tolerance of ", tolerance,"\n",omxRes@output$minimum," vs ",ppmlRes@output$minimum)	


	if (diff > tolerance)
	{
		print(sprintf("PPML found a lower minimum than pure numerical optimization. Diff %g", diff))
		return(TRUE)
	}

	print(sprintf("PPML equal to OpenMx minimum within tolerance. Diff %g", diff))

	return(FALSE)
#		omxCheckCloseEnough(res1@output$Minus2LogLikelihood, res2@output$Minus2LogLikelihood, 0.05)
}


# Functions to test PPML solutions against numerically optimized solutions
# TODO: Move to imxPPMLTester.R (?)
##' imxPPML.Test.Test
##'
##' Test that PPML solutions match non-PPML solutions.
##' 
##' @param model the MxModel to evaluate
##' @param checkLL whether to check log likelihood
##' @param checkByName check values using their names
##' @param tolerance closeness tolerance for check
##' @param testEstimates whether to test for the same parameter estimates
##' @details
##' This is an internal function used for comparing PPML and non-PPML solutions.
##' Generally, non-developers will not use this function.
imxPPML.Test.Test <- function(model, checkLL = TRUE, checkByName = FALSE, tolerance=0.5, testEstimates=TRUE) {
	# TODO: Wrap in timing functions for profiling
	model@options$UsePPML = "No"
	res1 <- mxRun(model, suppressWarnings = TRUE) # Standard fit
	model@options$UsePPML = "Yes"
	res2 <- mxRun(model, suppressWarnings = TRUE)	# PPML fit
	# /TODO

	# First model not transformed
	didNotUsePPMLon1 <- is.null(res1@options$UsePPML) || (res1@options$UsePPML == 'No')
	omxCheckTrue( didNotUsePPMLon1 )

	# Second model transformed
	usedPPMLon2 <- (!is.null(res2@options$UsePPML) &&
	 !(res2@options$UsePPML == "Inapplicable" || res2@options$UsePPML == "Yes" || res2@options$UsePPML == "No" ))
	omxCheckTrue( usedPPMLon2 )


	
	# NOTE: Not checking Hessians, they're never very close -- not sure that
	# a comparison is actually meaningful

	# NOTE: checkLL parameter can turn off log likelihood checking for this test,
	# useful for non-homogeneous error variance cases where the transformed log 
	# likelihood is different.
	if (testEstimates)
		PPML.Test.CheckFits(res1, res2, tolerance=tolerance, checkLL = checkLL, checkByName = checkByName, checkHessians = FALSE) # Check standard fit vs PPML model
	else
		PPML.Tool.CheckPPMLDidEqualOrBetter(res1, res2, tolerance)
	
}

PPML.Test.CheckFits <- function(res1, res2, tolerance, checkHessians = TRUE, checkLL, checkByName) {
	# Check -2logLLs versus each other
	# PPML must be <= OpenMx version, w/in tolerance
	if (checkLL)
	{
		PPMLDidBetter <- PPML.Tool.CheckPPMLDidEqualOrBetter(res1, res2, tolerance)	
		if (PPMLDidBetter)
			return() # Estimates, etc are not going to be the same if PPML found a different minimum
	}	
		
	if (checkByName) {
		
		for (label in labels(res1@output$estimate)) {
			omxCheckCloseEnough(res1@output$estimate[label], res2@output$estimate[label], tolerance)
		}
		
		oneInTwo <- lapply(labels(res1@output$estimate), function(l) { which(labels(res2@output$estimate) == l) })
		
		omxCheckCloseEnough(res1@output$standardErrors, as.matrix(res2@output$standardErrors[unlist(oneInTwo)]), tolerance)
		
	} else {
		# Check parameters versus each other
		# Parameters will be in the same order, so a simple iteration should work
		for ( i in 1:length(res1@output$estimate) ) {
			omxCheckCloseEnough(res1@output$estimate[i], res2@output$estimate[i], tolerance)
		}
		
		# Check standard errors versus each other
		# Should just be able to check the two vectors directly against each other
#		if (!single.na(res1@output$standardErrors) && !single.na(res2@output$standardErrors))
	#		omxCheckWithinPercentError(res1@output$standardErrors, res2@output$standardErrors, 1) # 1% difference allowed for
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

##' imxPPML.Test.Battery
##'
##' PPML can be applied to a number of special cases.  This function will test the given model for
##' all of these special cases.
##'
##' Requirements for model passed to this function:
##' - Path-specified
##' - Means vector must be present
##' - Covariance data (with data means vector)
##' - (Recommended) All error variances should be specified on the
##' diagonal of the S matrix, and not as a latent with a loading only
##' on to that manifest
##'
##' Function will test across all permutations of:
##' - Covariance vs Raw data
##' - Means vector present vs Means vector absent
##' - Path versus Matrix specification
##' - All orders of permutations of latents with manifests
##'
##' @param model the model to test
##' @param verbose whether to print diagnostics
##' @param testMissingness try with missingness
##' @param testPermutations try with permutations
##' @param testEstimates examine estimates
##' @param testFakeLatents try with fake latents
##' @param tolerances a vector of tolerances
imxPPML.Test.Battery <- function(model, verbose=FALSE, testMissingness = TRUE, testPermutations=TRUE, testEstimates=TRUE, testFakeLatents=TRUE, tolerances=c(.001, .001, .001))
{
	# Cov data check
	if (!model$data@type == "cov")
		return(NULL)
	# Data means vector presence check
	if (is.null(model$data@means))
		return(NULL)
	# Means vector presence check
	if (is.na(model$expectation@M))
		return(NULL)
	# Path-specified check
	if ( is.null(dimnames(model$S)) || !is.na(model$expectation@dims) )
		return(NULL)

	if (testMissingness)
		nM <- 1
	else
		nM <- 0

	

	# No missingness, missingness
	for ( missingness in 0:nM)
	{
		if (as.logical(missingness))
			if (verbose) print("Testing with missingness: ")

		# Test with no fake latents
		if (verbose) print("Testing with no fake latents: ")
		testModel <- model
		PPML.Test.Battery.LowLevel(testModel, verbose, missingness=as.logical(missingness), tolerances=tolerances, testPermutations=testPermutations, testEstimates=testEstimates)
		gc()

		# Test with varying numbers of fakeLatents
		if (testFakeLatents) {
			for (i in 1:length(model@manifestVars))
			{ 
				if (verbose) print(sprintf("Testing with %d fake latents: ", i))
				testModel@options$UsePPML <- "Yes" # Needs to be enabled for fakeLatentFoldout to get solveType
				testModel <- PPML.Tool.fakeLatentFoldout(testModel, i)
				PPML.Test.Battery.LowLevel(testModel, verbose=verbose, missingness=as.logical(missingness), tolerances=tolerances, testPermutations=testPermutations, testEstimates=testEstimates	)
				gc()
			}
		}		
	}
}

# This function does bitwise counting up to nMans bits
# Skips 0 == 0,0,...,0 pattern
PPML.Tool.EnumerateMissingnessPatterns <- function(nMans)
{
	patt <- logical(nMans)
	patt[1] <- TRUE
	
	M <- rbind(patt)
	while(TRUE)
	{
		for (i in 1:nMans)
		{
			if (!patt[i])
			{
				patt[i] <- TRUE

				if (i > 1)
					patt[1:i-1] <- FALSE

				M <- rbind(M, patt)
				break
			}
		}
		if (all(patt))
			break
	}
	return(M)
}


PPML.Test.Battery.LowLevel <- function(model, verbose = FALSE, missingness = FALSE, tolerances, testPermutations, testEstimates) {

	# Function only called by imxPPML.Test.Battery
	# --> Guaranteed to be path-specified
	
	# Get tolerances
	covMTolerance <- tolerances[1]
	covTolerance <- tolerances[2]
	rawTolerance <- tolerances[3]
	
	# Determine flags -- NA tolerance -> Don't check
	checkCovM <- !single.na(covMTolerance)
	checkCov <- !single.na(covTolerance)
	checkRaw <- !single.na(rawTolerance)

	# Create raw data
	# Pulled this out: means=as.vector(model$data@means),
	# Don't specify means for raw data
	rawData <- mxData(type="raw", numObs = model$data@numObs, mvrnorm(n=model$data@numObs, mu = model$data@means, model$data@observed) )
	
	if (missingness)
	{
		# Number of patterns = 2^num_mans - 1  (-1 because there is no data-less pattern)
		numPatts <- 2^length(model@manifestVars) - 1

		# Length of each missingness pattern in the dataset
		pattLen <- floor(model$data@numObs / numPatts)

		if (pattLen == 0)
			return(NULL) # Not enough observations to do all missingness patterns, abort

		patts <- PPML.Tool.EnumerateMissingnessPatterns(length(model@manifestVars))
		
		for ( i in 0:(dim(patts)[[1]]-1) )
		{
			patt <- as.vector(patts[i+1, ])
			
			for ( j in 1:length(patt) )
			{
				for ( k in 1:pattLen )
				{
					if (!patt[j])
						rawData@observed[i*pattLen + k, j] <- NA
				}
			}
		}
	
	}
	
	
	#### Path
	if (verbose) print("Testing Path-specified Models:")
	testModel <- model

	### +Expected Means
	## Cov
	if (checkCovM && !missingness)
	{
		if (verbose) print("Testing Covariance Data w/ Expected Means...")
		imxPPML.Test.Test(testModel, tolerance=covMTolerance, testEstimates=testEstimates)
	}
	
	## Raw
	if (checkRaw)
	{
		if (verbose) print("Testing Raw Data w/ Expected Means...")
		imxPPML.Test.Test(mxModel(testModel, rawData), tolerance=rawTolerance, testEstimates=testEstimates)
	}
		
	### -Expected Means
	if (checkCov && !missingness)
	{
		testModel[[testModel$expectation@M]] <- NULL
		testModel$expectation@M <- as.character(NA)
		## Cov
		if (verbose) print("Testing Covariance Data w/o Expected Means...")
		imxPPML.Test.Test(mxModel(testModel, mxData(type=testModel$data@type, numObs=testModel$data@numObs, observed=testModel$data@observed)), tolerance=covTolerance, testEstimates=testEstimates)
		## Raw: Can't test for Raw -Expected Means, as raw data requires an expected means vector
	}
	
	
	#### Matrix
	if (verbose) print("Testing Matrix-specified Models:")
	
	matrixModel <- PPML.Tool.PathToMatrix(model)
	# Try each unique permutation

	# -First permutation in the array is just the unpermuted case
	# -OpenMx doesn't like having its matrices permuted, so just get the unPPMLed
	# output from the unpermuted case, and permute it accordingly to check
	# against PPMLed results
	
	testModel <- matrixModel
	if (verbose) print("Testing unpermuted model...")
	### +Means
	## Cov
	if (checkCovM && !missingness)
	{
		if (verbose) print("Testing Covariance Data w/ Expected Means...")
		testModel@options$UsePPML <- "No"
		vanillaCovMRes <- mxRun(testModel)
		testModel@options$UsePPML <- "Yes"
		PPMLCovMRes <- mxRun(testModel)

		usedPPML <- (!is.null(PPMLCovMRes@options$UsePPML) &&
		 !(PPMLCovMRes@options$UsePPML == "Inapplicable" || PPMLCovMRes@options$UsePPML == "Yes" || PPMLCovMRes@options$UsePPML == "No" ))
		omxCheckTrue( usedPPML )
		didBetter <- PPML.Tool.CheckPPMLDidEqualOrBetter(vanillaCovMRes, PPMLCovMRes, tolerance=covMTolerance)
		if (!didBetter && testEstimates)
			omxCheckCloseEnough(vanillaCovMRes@output$estimate, PPMLCovMRes@output$estimate, covMTolerance)
	}


	## Raw
	if (checkRaw)
	{
		if (verbose) print("Testing Raw Data w/ Expected Means...")
		testModel@options$UsePPML <- "No"
		vanillaRawRes <- mxRun(mxModel(testModel, rawData) )
	
	
		testModel@options$UsePPML <- "Yes"
		PPMLRawRes <- mxRun(mxModel(testModel, rawData))
		usedPPML <- (!is.null(PPMLRawRes@options$UsePPML) &&
		 !(PPMLRawRes@options$UsePPML == "Inapplicable" || PPMLRawRes@options$UsePPML == "Yes" || PPMLRawRes@options$UsePPML == "No" ))
		omxCheckTrue( usedPPML )
		didBetter <- PPML.Tool.CheckPPMLDidEqualOrBetter(vanillaRawRes, PPMLRawRes, tolerance=rawTolerance)
		if (!didBetter && testEstimates)
			omxCheckCloseEnough(vanillaRawRes@output$estimate, PPMLRawRes@output$estimate, rawTolerance)	
	}

	### -Means
	if (checkCov && !missingness)
	{
		testModel[[testModel$expectation@M]] <- NULL
		testModel$expectation@M <- as.character(NA)
		dataNoMeans <- mxData(type=testModel$data@type, numObs=testModel$data@numObs, observed=testModel$data@observed)
		## Cov
		if (verbose) print("Testing Covariance Data w/o Expected Means...")
		testModel@options$UsePPML <- "No"
		vanillaCovRes <- mxRun(mxModel(testModel, data=dataNoMeans))

		testModel@options$UsePPML <- "Yes"
		PPMLCovRes <- mxRun(mxModel(testModel, data=dataNoMeans))
		usedPPML <- (!is.null(PPMLCovRes@options$UsePPML) &&
		 !(PPMLCovRes@options$UsePPML == "Inapplicable" || PPMLCovRes@options$UsePPML == "Yes" || PPMLCovRes@options$UsePPML == "No" ))
		omxCheckTrue( usedPPML )
		didBetter <- PPML.Tool.CheckPPMLDidEqualOrBetter(vanillaCovRes, PPMLCovRes, tolerance=covTolerance)
		if (!didBetter && testEstimates)
			omxCheckCloseEnough(vanillaCovRes@output$estimate, PPMLCovRes@output$estimate, covTolerance)
	}
	

	if (!testPermutations)
		return();

	# Want to check all relevant permutations of manifests and latents; however,
	# switching the order of any two latents or any two manifests should not pose any
	# kind of problem.  So, do only unique permutations treating all manifests as 
	# the same and all latents as the same for the purposes of the permutation.
	defPerm <- 1:dim(model[[model$expectation@F]])[[2]]
	manOrLat <- apply(model[[model$expectation@F]]@values, 2, sum)
	manIndices <- which(manOrLat == 1)
	latIndices <- which(manOrLat == 0)
	# Get unique permutations
	perms <- unique(PPML.Tool.EnumeratePermutations(manOrLat))
	# Insert actual indices for mans or lats in to permutations
	perms <- t(apply(perms, 1, function(row) { row[which(row==1)] <- manIndices; row[which(row==0)] <- latIndices; return(row) } ))


	# PERMUTATIONS OF MATRIX MODELS
	# Get a standard ordering for estimates
#	estNamesM <- names(vanillaCovMRes@output$estimate)
#	estNames <- names(vanillaCovRes@output$estimate)

	# Minima tracking
	covMEsts <- list()
	covEsts <- list()
	rawEsts <- list()

	for (i in 2:dim(perms)[[1]]) {

		if (verbose)
		{
			print("Testing permutation: ")
			print(as.vector(perms[i,]))
		} 
		testModel <- matrixModel
		# PERMUTING THE MODEL
		# Build permutation matrix
		#P <- matrix(nrow=length(defPerm), ncol=length(defPerm),
		#	unlist(lapply(perms[i,], function(item){ as.numeric(defPerm == item) } )) )
		
		# Permute matrices
		# A %*% P
		testModel[[testModel$expectation@A]] <- testModel[[testModel$expectation@A]][ , perms[i, ] ]
		# t(P) %*% (A %*% P)
		testModel[[testModel$expectation@A]] <- testModel[[testModel$expectation@A]][ perms[i, ], ]
		
		# For S -- Special considerations necessary, as the operating on S in this manner breaks
		# the symmetry of the matrix and changes the output
		#   So, instead of directly: Extract labels, free, values, operate on them, reinsert
		Sfree <- testModel[[testModel$expectation@S]]@free
		Slabels <- testModel[[testModel$expectation@S]]@labels
		Svalues <- testModel[[testModel$expectation@S]]@values
		# S %*% P
		Sfree <- Sfree[ , perms[i, ] ]
		Slabels <- Slabels[ , perms[i, ] ]
		Svalues <- Svalues[ , perms[i, ] ]
		# t(P) %*% (S %*% P)
		Sfree <- Sfree[ perms[i, ], ]
		Slabels <- Slabels[ perms[i, ], ]
		Svalues <- Svalues[ perms[i, ], ]
		# Reinsert
		testModel[[testModel$expectation@S]]@free <- Sfree
		testModel[[testModel$expectation@S]]@labels <- Slabels
		testModel[[testModel$expectation@S]]@values <- Svalues
		
		# F %*% P
		testModel[[testModel$expectation@F]] <- testModel[[testModel$expectation@F]][ , perms[i, ] ]

		# M %*% P
		testModel[[testModel$expectation@M]] <- testModel[[testModel$expectation@M]][ , perms[i, ] ]


		# Don't need to permute data because manifest variables are always
		# in the same order.
		# Permute data cov matrix, means vector
#		dcov <- testModel$data@observed
#		dcov <- dcov[ , intersect(perms[i, ], manIndices) ]
#		dcov <- dcov[ intersect(perms[i, ], manIndices) , ]
#		testModel$data@observed <- dcov
#		dmeans <- t(as.matrix(testModel$data@means[ intersect(perms[i, ], manIndices) ]))
#		dimnames(dmeans)[[2]] <- (dimnames(testModel$data@means)[[2]])[ intersect(perms[i,], manIndices) ]
#		testModel$data@means <- dmeans
		

		# Permute raw data
#		rawDataPerm <- rawData
#		rawDataPerm@observed <- rawDataPerm@observed[ , intersect(perms[i, ], manIndices) ]

		# Permute expectation dims
		testModel$expectation@dims <- testModel$expectation@dims[ perms[i, ] ]
		
		# Enable PPML for testing permuted model
		testModel@options$UsePPML <- "Yes"

		# NOTE: Estimates are not checked for the permutated models
		# Permuting the matrices in this way has a way of finding alternative, equally good minimums
		# Should probably investigate this tendency further

		### +Means
		## Cov
		if (checkCovM && !missingness)
		{
			if (verbose) print("Testing Covariance Data w/ Expected Means...")
			PPMLCovMRes <- mxRun(testModel)
			usedPPML <- (!is.null(PPMLCovMRes@options$UsePPML) &&
			 !(PPMLCovMRes@options$UsePPML == "Inapplicable" || PPMLCovMRes@options$UsePPML == "Yes" || PPMLCovMRes@options$UsePPML == "No" ))
			omxCheckTrue( usedPPML )
			didBetter <- PPML.Tool.CheckPPMLDidEqualOrBetter(vanillaCovMRes, PPMLCovMRes, covMTolerance)
			#omxCheckCloseEnough(vanillaCovMRes@output$estimate, PPMLCovMRes@output$estimate, .05)
#			covMEsts <- append(covMEsts, list(PPMLCovMRes@output$estimate[estNamesM]))
		}

		## Raw
		if (checkRaw)
		{
			if (verbose) print("Testing Raw Data w/ Expected Means...")
			PPMLRawRes <- mxRun(mxModel(testModel, rawData))
			usedPPML <- (!is.null(PPMLRawRes@options$UsePPML) &&
			 !(PPMLRawRes@options$UsePPML == "Inapplicable" || PPMLRawRes@options$UsePPML == "Yes" || PPMLRawRes@options$UsePPML == "No" ))
			omxCheckTrue( usedPPML )
			didBetter <- PPML.Tool.CheckPPMLDidEqualOrBetter(vanillaRawRes, PPMLRawRes, rawTolerance)
			#omxCheckCloseEnough(vanillaRawRes@output$estimate, PPMLRawRes@output$estimate, .05)	
			#rawEsts <- append(rawEsts, list(PPMLRawRes@output$estimate[estNamesM]))
		}
		
		### -Means
		if (checkCov && !missingness)
		{
			testModel[[testModel$expectation@M]] <- NULL
			testModel$expectation@M <- as.character(NA)
			## Cov
			if (verbose) print("Testing Covariance Data w/o Expected Means...")
			PPMLCovRes <- mxRun(mxModel(testModel, data=dataNoMeans))
			usedPPML <- (!is.null(PPMLCovRes@options$UsePPML) &&
			 !(PPMLCovRes@options$UsePPML == "Inapplicable" || PPMLCovRes@options$UsePPML == "Yes" || PPMLCovRes@options$UsePPML == "No" ))
			didBetter <- PPML.Tool.CheckPPMLDidEqualOrBetter(vanillaCovRes, PPMLCovRes, covTolerance)
#			covEsts <- append(covEsts, list(PPMLCovRes@output$estimate[estNames]))
			#omxCheckCloseEnough(vanillaCovRes@output$estimate, PPMLCovRes@output$estimate, .05)
		}
	}
}

PPML.Tool.EnumeratePermutations <- function(elements) {
	if (length(elements) == 1) 
		return(elements)
	pMat <- matrix(nrow=factorial(length(elements)), ncol=length(elements))
	blockSize <- factorial(length(elements) - 1)
	for (i in 0:(length(elements)-1) ) {
		pMat[ (i*blockSize+1) : ((i+1)*blockSize), 1] <- elements[i+1]
		pMat[ (i*blockSize+1) : ((i+1)*blockSize), 2:length(elements)] <- as.matrix(PPML.Tool.EnumeratePermutations(elements[-(i+1)]))
	}
	return(pMat)
}

# Function to respecify a path-specified model as a matrix-specified model
PPML.Tool.PathToMatrix <- function(model) {
	### Extract expectation
	expectation <- model$expectation
	if(is.null(expectation) || !is(expectation, "MxExpectationRAM")) {
		return(NA)
	}
	
	### Safely extract ASF(M) matrices by using names from the objective
	Aname <- expectation@A
	Sname <- expectation@S
	Fname <- expectation@F
	Mname <- expectation@M
	
	Amatrix <- model[[Aname]]
	Smatrix <- model[[Sname]]
	Fmatrix <- model[[Fname]]
	
	# Extract variable names for new expectation
	varNames <- dimnames(Fmatrix)[[2]]
	
	### Remove dimnames from matrices
	dimnames(Amatrix) <- NULL
	dimnames(Smatrix) <- NULL
	dimnames(Fmatrix) <- NULL
	
	# Special: Matrix of expected means
	Mmatrix <- NULL
	if (!single.na(Mname)) {
		Mmatrix <- model[[Mname]]
		dimnames(Mmatrix) <- NULL
	}

	# Create new expectation
	newexpectation <- mxExpectationRAM(A=Aname,S=Sname,F=Fname,M=Mname, dimnames=varNames)
	# Remove manifestvars, latentvars
	model@latentVars <- character(0)
	model@manifestVars <- character(0)

	# Return matrix-specified model
	return( mxModel(name = model@name, Amatrix, Smatrix, Fmatrix, Mmatrix, newexpectation, model$fitfunction, model$data, model@constraints ) )
}

# This function only accepts path-specified models
# Takes the manifest variable at manIndex, and folds its variance out in to a fake latent
# Fake latent := Manifest variance specified by a latent variable with only one loading
#  of value 1, to the manifest in question ==> Variance of the manifest is specified
#  in the fakelatent, not in the manifest variable itself 
PPML.Tool.fakeLatentFoldout <- function(model, manIndex)
{
	# Path-specified check
	if ( is.null(dimnames(model$S)) || !single.na(model$expectation@dims) )
		return(NULL)

	Sname <- model$expectation@S
	Smatrix <- model[[Sname]]

	
	# Idiot check: Make sure manifest[manIndex] exists
	if ( manIndex > length(model@manifestVars) )
		return(NULL)

	man <- model@manifestVars[manIndex]

	# Save info
	FLvar <- Smatrix@values[man,man]
	FLlabel <- Smatrix@labels[man,man]
	FLfree <- Smatrix@free[man,man]

	# Remove normally defined variance from Smatrix
	Smatrix@values[man,man] <- 0
	Smatrix@free[man,man] <- FALSE
	Smatrix@labels[man,man] <- NA

	# Reinsert stripped Smatrix
	model[[Sname]] <- Smatrix	

	## Add fakeLatent
	# Generate name -- man is the name of the variable in path-specified models
	FLname <- sprintf("%s_FL", man)
	
	# Create path from fakelatent to manifest variable
	oneHeadedPath <- mxPath(from=FLname, to=man, arrows=1, free=FALSE, values=1)
	
	# Create two-headed arrow on fake latent
	twoHeadedPath <- mxPath(from=FLname, arrows=2, free=FLfree, values = FLvar, labels = FLlabel)

	# Recreate model
	model <- mxModel(model, latentVars=FLname, oneHeadedPath, twoHeadedPath)
	return(model)	
}
