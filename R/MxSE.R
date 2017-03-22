#   Copyright 2007-2017 The OpenMx Project
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
##' @details
##' x can be the name of an algebra, a bracket address, named entity
##' or arbitrary expression. It is a frontend-only file that works
##' much like mxEval. When the \code{details} argument is TRUE, the full
##' sampling covariance matrix of \code{x} is also returned as part of a list.
##' The square root of the diagonals of this sampling covariance matrix are
##' the standard errors.
##'
##' @param x the parameter to get SEs on (reference or expression)
##' @param model the \code{\link{mxModel}} to use.
##' @param details logical. Whether to provide further details, e.g. the full sampling covariance matrix of x.
##' @param ... further named arguments passed to \code{\link{mxEval}}
##' @return SE value(s) returned as a matrix when \code{details} is FALSE.
##' When \code{details} is TRUE, a list of the SE value(s) and the full sampling covariance matrix.
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
##' mxSE('lambda5', model = m1)
##' mxSE(lambda1^2, model = m1)
mxSE <- function(x, model, details=FALSE, ...){
	isCallEtc <- any(c('call', 'language', 'MxAlgebraFormula') %in% is(match.call()$x))
	if(isCallEtc){
		message('Treating first argument as an expression')
		xalg <- mxAlgebraFromString(Reduce(paste, deparse(match.call()$x)), name='onTheFlyAlgebra')
		x <- "onTheFlyAlgebra"
		model <- mxModel(model, xalg)
	} else if('character' %in% is(x)){
		message('Treating first argument as character named entity in the model')
	} else {
		stop("Please, sir.  'x' must be either the name of an entity in the model, or an expression for an MxAlgebra.")
	}
	
	sefun <- function(x = NULL, model, alg){
		if(is.null(x)){x <- omxGetParameters(model)}
		param <- omxGetParameters(model)
		paramNames <- names(param)
		model <- omxSetParameters(model, values=x, labels=paramNames, free=TRUE)
		return(c(mxEvalByName(alg, model, compute=TRUE)))
	}
	
	# Get current algebra/matrix values:
	freeparams <- omxGetParameters(model)
	paramnames <- names(freeparams)
	matrix(NA)
	zoutMat <- try(mxEvalByName(x, model, compute=TRUE))
	if(class(zoutMat) %in% "try-error"){
		stop(paste0("I had a problem evaluating your expression. Check that mxEval works for it.\nRecall that elements of submodels are addressed as submodelName.objectName\nFor example, 'sub1.bob' refers to the object 'bob' in the submodel named 'sub1'."))
	}
	zoutVec <- sefun(x=freeparams, model=model, alg=x)
	
	if(length(model@output) > 0){
		if(model@output$infoDefinite){
			# solve() will fail if Hessian is computationally singular;
			# chol2inv() will still fail if Hessian is exactly singular.
			ParamsCov <- 2*chol2inv(chol(model@output$hessian))
			dimnames(ParamsCov) <- dimnames(model@output$hessian)
		} else{
			# An indefinite Hessian usually means some SEs will be NaN:
			ParamsCov <- 2*solve(model@output$hessian)
		}
	} else {
		stop("Model does not have output.  I'm a doctor, not a bricklayer!\nWas this model run with mxRun?")
	}
	
	covParam <- ParamsCov[paramnames,paramnames] # <-- submodel will usually not contain all free param.s
	jacTrans <- numDeriv::jacobian(func = sefun, x = freeparams, model = model, alg = x)
	covSparam <- jacTrans %*% covParam %*% t(jacTrans)
	# dimnames(covSparam) <- list(rownames(zoutMat), colnames(zoutMat))
	SEs <- sqrt(diag(covSparam))
	SEsMat <- matrix(SEs, nrow = nrow(zoutMat), ncol = ncol(zoutMat))
	if(details==TRUE){
		return(list(SE=SEsMat, Cov=covSparam))
	} else{
		return(SEsMat)
	}
}

