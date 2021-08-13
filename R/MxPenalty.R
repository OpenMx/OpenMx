#
#   Copyright 2021 by the individuals mentioned in the source code history
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

#' MxPenalty
#'
#' This is an internal class and should not be used directly.
#'
#' @aliases
#' $,MxPenalty-method
#' $<-,MxPenalty-method
#' [[<-,MxPenalty-method
setClass(Class = "MxPenalty",
          representation = representation(
           name = "character",
           type = "character",
           params = "MxCharOrNumber",
           epsilon = "numeric",
           scale = "numeric",
           smoothProportion = "numeric",
           hyperparameters = "MxCharOrNumber",
           hpranges = "list"
          )
         )

setMethod("initialize", "MxPenalty",
          function(.Object, name, type, reg_params, epsilon, scale,
                   smoothProportion, hyperparams=NULL, hpranges=NULL) {
            .Object@name <- name
            .Object@type <- type
            .Object@params <- reg_params
            .Object@epsilon <- epsilon
            .Object@scale <- scale
            .Object@smoothProportion <- smoothProportion
            .Object@hyperparameters <- hyperparams
            .Object@hpranges <- hpranges
            return(.Object)
          }
)

setMethod("$", "MxPenalty", imxExtractSlot)

setReplaceMethod("$", "MxPenalty",
	function(x, name, value) {
		return(imxReplaceSlot(x, name, value, check=TRUE))
	}
)

setMethod("names", "MxPenalty", slotNames)

##' imxRegularizationTypes
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @details Types of regularization penalties.
imxRegularizationTypes <- c('lasso', 'ridge', 'elasticNet')

##' This function creates a regularization penalty object
##'
##' @param how what kind of regularization function to use
##' @param hyperparams a character vecot of hyperparameter names
##' @param hpranges a named list of hyperparameter ranges. Used in search if no ranges are specified.
##' @template args-regularize
##' @param smoothProportion what proportion of the region between \code{epsilon} and zero should be used to smooth the penalty function
##'
##' @details \code{mxRegularize} expects to find an \link{mxMatrix} with
##'   free parameters that correspond to all named hyperparameters.
##'
##' Gradient descent optimizers are designed for and work best on
##' smooth functions.  All of the regularization penalties implemented
##' traditionally contain discontinuities.  By default, OpenMx uses
##' smoothed versions of these functions. Smoothing is controlled by
##' \code{smoothProportion}. If \code{smoothProportion} is zero then
##' the traditional discontinuous functions are used. Otherwise,
##' \code{smoothProportion} of the region between \code{epsilon} and
##' zero is used for smoothing.
mxRegularize <- function(what, epsilon=1e-5, scale=1,
                         how=imxRegularizationTypes,
                         smoothProportion = 0.05,
                         hyperparams=c(),
                         hpranges=list(),
                         name=NULL) {
    how <- match.arg(how)
    if(is.null(name)) name <- imxUntitledName()
    if(length(grep("data[\\.]", what))>0) {
      stop("Error in regularization penalty ", name,
           ": cannot regularize definition vars")
    }
    if (length(what) %% length(epsilon) != 0) {
      stop("length(what) %% length(epsilon) != 0")
    }
    if (length(what) %% length(scale) != 0) {
      stop("length(what) %% length(scale) != 0")
    }
    new("MxPenalty", name=name, type=how, reg_params=what, epsilon=epsilon, scale=scale,
        smoothProportion=smoothProportion, hyperparams=hyperparams, hpranges=hpranges)
}

##' MxRegularizeLASSO
##'
##' Least Absolute Selection and Shrinkage Operator regularization
##'
##' @param lambda strength of the penalty to be applied at starting values (default 0)
##' @param lambda.step step function for lambda step (default .01)
##' @param lambda.max end of lambda range (default .4)
##' @param lambda.min minimum lambda value (default lambda)
##' @template args-regularize
##'
#' @export
mxRegularizeLASSO <- function(what, name, lambda=0, lambda.step=.01, lambda.max=NA, lambda.min=NA,
                              epsilon=1e-5, scale=1) {
  if(is.na(lambda.min)) lambda.min=lambda
  if(is.na(lambda.max)) lambda.max=lambda.min+40*lambda.step
  mxRegularize(what, epsilon, scale, how="lasso", hyperparams=c('lambda'), hpranges=list(lambda=seq(lambda, lambda.max, by=lambda.step)), name=name)
}

##' MxRegularizeRidge
##'
##' Ridge regression regularization
##'
##' @param lambda strength of the penalty to be applied at start (default 0)
##' @param lambda.step lambda step during penalty search (default 0.01)
##' @param lambda.max when to end the lambda search (default 0.4)
##' @param lambda.min minimum lambda value (default lambda)
##' @template args-regularize
##'
#' @export
mxRegularizeRidge <- function(what, name, lambda=0, lambda.step=.01, lambda.max=.4, lambda.min=NA,
                              epsilon=1e-5, scale=1) {
  if(is.na(lambda.min)) lambda.min=lambda
  mxRegularize(what, epsilon, scale, how="ridge", hyperparams=c('lambda'),
               hpranges = list(lambda=seq(lambda, lambda.max, lambda.step)), name=name)
}

##' MxRegularizeElasticNet
##'
##' Elastic net regularization
##'
##' @param alpha strength of the mixing parameter to be applied at start (default 0.5).  Note that 0 indicates a ridge regression with penalty \deqn{\frac{lambda}{2}}{lambda / 2}, and 1 indicates a LASSO regression with penalty lambda.
##' @param alpha.step alpha step during penalty search (default 0.1)
##' @param alpha.max when to end the alpha search (default 1)
##' @param lambda strength of the penalty to be applied at starting values (default 0)
##' @param lambda.step step function for lambda step (default .01)
##' @param lambda.max end of lambda range (default .4)
##' @param lambda.min beginning of the lambda range (default lambda)
##' @param alpha.min beginning of the alpha range (default 0)
##' @template args-regularize
##'
##' @details Applies elastic net regularization.  Elastic net is a weighted combination of ridge and LASSO penalties.
##'
#' @export
mxRegularizeElasticNet <- function(what, name,
                                   alpha=0,  alpha.step=.1,  alpha.max=1,
                                   lambda=0, lambda.step=.1, lambda.max=.4,
                                   alpha.min=NA, lambda.min=NA,
                                   epsilon=1e-5, scale=1) {
  if(is.na(lambda.min)) lambda.min=lambda
  if(is.na(alpha.min)) alpha.min=alpha
  mxRegularize(what, epsilon, scale, how="elasticNet", hyperparams=c('alpha', 'lambda'),
               hpranges = list(alpha=seq(alpha, alpha.max, alpha.step), lambda=seq(lambda, lambda.max, by=lambda.step)), name=name)
}

generatePenaltyList <- function(flatModel, modelname, parameters, labelsData) {
	got <- lapply(flatModel@regularizations, function(p1) {
    hpnames <- p1@hyperparameters
    p1@hyperparameters <- match(p1@hyperparameters, names(parameters)) - 1L
    if (any(is.na(p1@hyperparameters))) {
      stop(paste("Cannot locate regularization parameters",
                 omxQuotes(hpnames[is.na(p1@hyperparameters)])))
    }
    pnames <- p1@params
    p1@params <- match(p1@params, names(parameters)) - 1L
    if (any(is.na(p1@params))) {
      # maybe too strict?
      stop(paste("Cannot locate parameters",
                 omxQuotes(pnames[is.na(p1@params)]), "to regularize"))
    }
    p1
  })
  if (length(got)) {
    totalParams <- sum(sapply(flatModel@regularizations, function(p1) length(p1@params)))
    if (totalParams != length(Reduce(union, lapply(flatModel@regularizations, function(p1) p1@params))))
      stop("Attempt to regularize the same parameter more than once")
  }
  got
}
