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
# Date: 2016-08-12 09:28:02 by mhunter
# Date: 2016-08-13 add roxygen by tbates
# Filename: MxSE.R
# Purpose: Define the mxSE() function.		
#  This function lets the user obtain standard errors for arbitrary		
#  expressions, named entities, and algebras.		
#  It is a frontend-only file that works much like mxEval.
#------------------------------------------------------------------------------


##' Compute standard errors in OpenMx
##' 
##' @description
##' This function allows you to obtain standard errors for arbitrary
##' expressions, named entities, and algebras.
##' 
##' @param x the parameter to get SEs on (reference or expression)
##' @param model the \code{\link{mxModel}} to use.
##' @param details logical. Whether to provide further details, e.g. the full
##' sampling covariance matrix of x.
##' @param cov optional matrix of covariances among the free parameters. If 
##' missing, the inverse Hessian from the fitted model is used.
##' @param forceName logical; defaults to \code{FALSE}.  Set to \code{TRUE}
##' if \code{x} is an R symbol that refers to a character string.
##' @param silent logical; defaults to \code{FALSE}.  If \code{TRUE},
##' message-printing is suppressed.
##' @param ... further named arguments passed to \code{\link{mxEval}}
##' @param defvar.row which row to load for any definition variables
##' @param data name of data from which to load definition variables
##' 
##' @details
##' x can be the name of an algebra, a bracket address, named entity
##' or arbitrary expression.
##' When the \code{details} argument is TRUE, the full
##' sampling covariance matrix of \code{x} is also returned as part of a list.
##' The square root of the diagonals of this sampling covariance matrix are
##' the standard errors.
##' 
##' When supplying the \code{cov} argument, take care that the free parameter
##' covariance matrix is given, not the information matrix.  These 
##' two are inverses of one another.
##' 
##' This function uses the delta method to compute the standard error of arbitrary
##' and possibly nonlinear functions of the free parameters.  The delta method
##' makes a first-order Taylor approximation of the nonlinear function.  The
##' nonlinear function is a map from all the free parameters to some transformed
##' subset of parameters: the linearization of this map is given by the Jacobian
##' \eqn{J}.  In equation form, the delta method computes standard errors by the following:
##' 
##' \deqn{J^T C J}
##' 
##' where \eqn{J} is the Jacobian of the nonlinear parameter transformation
##' and \eqn{C} is the covariance matrix of the free parameters (e.g., two
##' times the inverse of the Hessian of the minus two log likelihood function).
##' 
##' @return SE value(s) returned as a matrix when \code{details} is FALSE.
##' When \code{details} is TRUE, a list of the SE value(s) and the full 
##' sampling covariance matrix.
##' @seealso - \code{\link{mxCI}}
##' @references - \url{https://en.wikipedia.org/wiki/Standard_error}
##' @examples
##' library(OpenMx)
##' data(demoOneFactor)
##' # ===============================
##' # = Make and run a 1-factor CFA =
##' # ===============================
##' 
##' latents  = c("G") # the latent factor
##' manifests = names(demoOneFactor) # manifest variables to be modeled
##' # ===========================
##' # = Make and run the model! =
##' # ===========================
##' m1 <- mxModel("One Factor", type = "RAM", 
##' 	manifestVars = manifests, latentVars = latents, 
##' 	mxPath(from = latents, to = manifests, labels=paste0('lambda', 1:5)),
##' 	mxPath(from = manifests, arrows = 2),
##' 	mxPath(from = latents, arrows = 2, free = FALSE, values = 1),
##' 	mxData(cov(demoOneFactor), type = "cov", numObs = 500)
##' )
##' m1 = mxRun(m1)
##' mxSE(lambda5, model = m1)
##' mxSE(lambda1^2, model = m1)
mxSE <- function(x, model, details=FALSE, cov, forceName=FALSE, silent=FALSE, ...,
		 defvar.row=as.integer(NA), data='data'){
  warnModelCreatedByOldVersion(model)
	if(length(model@output) > 0 && missing(cov)){
		ParamsCov <- try(vcov(model))
		if(is(ParamsCov,"try-error")){
			msg <- "Model does not have a reasonable vcov matrix or standard errors."
			stop(msg)
		}
		# if(length(model@output$infoDefinite) && !single.na(model@output$infoDefinite)){
		# 	# An indefinite Hessian usually means some SEs will be NaN:
		# 	ParamsCov <- 2*solve(model@output$hessian)
		# 	dimnames(ParamsCov) <- dimnames(model@output$hessian)
		# } else {
	} else if (missing(cov)){
		stop("Model does not have output and 'cov' argument is missing.  I'm a doctor, not a bricklayer!\nWas this model run with mxRun?")
	} else {
		ParamsCov <- cov
		if(is.null(dimnames(ParamsCov))){
			if(length(paramnames) == nrow(ParamsCov)){
				dimnames(ParamsCov) <- list(paramnames, paramnames)
			}else{
				stop(paste0("dimnames of user-supplied parameter covariance matrix are null\nand the number of rows (",  nrow(ParamsCov), ") do not match the number of free parameters (", length(paramnames), ")."))
			}
		}
	}
	
	xorig <- "x" #<--Initialize as something that will always be understandable in an error message.
	isCallEtc <- any(c('call', 'language', 'MxAlgebraFormula') %in% is(match.call()$x))
	ex <- try(eval(x), silent=TRUE)
	isChar <- !('try-error' %in% is(ex)) && is.character(ex)
	if(isCallEtc && !forceName && !isChar){
		if(!silent){message('Treating first argument as an expression')}
		xalg <- mxAlgebraFromString(Reduce(paste, deparse(match.call()$x)), name='onTheFlyAlgebra')
		xorig <- Reduce(paste, deparse(match.call()$x))
		x <- "onTheFlyAlgebra"
		model <- mxModel(model, xalg)
	} else if ('character' %in% is(x) && !isCallEtc) {
		if(!silent){message('Treating first argument as a character')}
		xalg <- mxAlgebraFromString(Reduce(paste, match.call()$x), name='onTheFlyAlgebra')
		xorig <- x
		x <- "onTheFlyAlgebra"
		model <- mxModel(model, xalg)
	} else if(isChar){
		if(!silent){message('Treating first argument as an object that stores a character')}
		xalg <- mxAlgebraFromString(ex, name='onTheFlyAlgebra')
		xorig <- ex
		x <- "onTheFlyAlgebra"
		model <- mxModel(model, xalg)
	} else {
		stop("Please, sir.  'x' must be either the name of an entity in the model, or an expression for an MxAlgebra.")
	}
	
	# Get current algebra/matrix values:
	freeparams <- omxGetParameters(model)
	paramnames <- names(freeparams)
	zoutMat <- try(mxEvalByName(x, model, compute=TRUE),silent=silent)
	if(is(zoutMat, "try-error")) {
		stop(paste0("Couldn't evaluate expression ", omxQuotes(xorig), ". Might help to check if it works in mxEval.\n",
		"Recall also that elements of submodels are addressed as submodelName.objectName\n",
		"For example, to refer to an object called 'bob' in submodel 'sub1', you would say 'sub1.bob'."))
	}
	
	covParam <- ParamsCov
	jModel <- mxModel(model, mxComputeJacobian(of=x, defvar.row=defvar.row, data=data))
	jModel <- mxRun(jModel, silent=TRUE)
	jacTrans <- jModel$compute$output$jacobian
	covSparam <- jacTrans %*% covParam %*% t(jacTrans)
	# dimnames(covSparam) <- list(rownames(zoutMat), colnames(zoutMat))
	if(any(diag(covSparam) < 0) || any(is.na(diag(covSparam)))){
		warning("Some diagonal elements of the repeated-sampling covariance matrix of the point estimates are less than zero or NA.\nI know, right? Set details=TRUE and check the 'Cov' element of this object.")
	}
	SEs <- suppressWarnings(sqrt(diag(covSparam)))
	SEsMat <- matrix(SEs, nrow = nrow(zoutMat), ncol = ncol(zoutMat))
	if(details==TRUE){
		return(list(SE=SEsMat, Cov=covSparam))
	} else{
		return(SEsMat)
	}
}

