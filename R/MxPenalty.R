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
           type = "character",
           params = "MxCharOrNumber",
           epsilon = "numeric",
           scale = "numeric",
           smoothProportion = "numeric",
           hyperparameters = "MxCharOrNumber",
           hpranges = "list",
           result = "matrix"
          ),
         contains = "MxBaseNamed")

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
            .Object
          }
)

setMethod("$", "MxPenalty", imxExtractSlot)

setReplaceMethod("$", "MxPenalty",
	function(x, name, value) {
        if(name == "result") {
            stop("You cannot set the result of a penalty.  Use mxRun() to populate the result, or mxEval() to compute it.")
        }
		return(imxReplaceSlot(x, name, value, check=TRUE))
	}
)

setMethod("names", "MxPenalty", slotNames)

##' imxPenaltyTypes
##'
##' This is an internal function exported for those people who know
##' what they are doing.
##'
##' @details Types of regularization penalties.
imxPenaltyTypes <- c('lasso', 'ridge', 'elasticNet')

##' This function creates a penalty object
##'
##' @param how what kind of function to use
##' @param hyperparams a character vector of hyperparameter names
##' @param hpranges a named list of hyperparameter ranges. Used in search if no ranges are specified.
##' @template args-regularize
##' @param smoothProportion what proportion of the region between \code{epsilon} and zero should be used to smooth the penalty function
##'
##' @details \code{mxPenalty} expects to find an \link{mxMatrix} with
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
mxPenalty <- function(what, epsilon=1e-5, scale=1,
                         how=imxPenaltyTypes,
                         smoothProportion = 0.05,
                         hyperparams=c(),
                         hpranges=list(),
                         name=NULL) {
    how <- match.arg(how)
    if(is.null(name)) name <- imxUntitledName()
    if(length(grep("data[\\.]", what))>0) {
      stop("Error in penalty ", name,
           ": cannot apply penalty to definition vars")
    }
    if (length(what) %% length(epsilon) != 0) {
      stop("length(what) %% length(epsilon) != 0")
    }
    if (length(what) %% length(scale) != 0) {
      stop("length(what) %% length(scale) != 0")
    }
    if (any(duplicated(what))) {
      stop(paste0("Parameters", omxQuotes(what[duplicated(what)]),
                  "added to penalty more than once"))
    }
    new("MxPenalty", name=name, type=how, reg_params=what, epsilon=epsilon, scale=scale,
        smoothProportion=smoothProportion, hyperparams=hyperparams, hpranges=hpranges)
}

regularizedToZeroParameters <- function(pen, model) {
  mask <- coef(model)[pen@params] / pen@scale <= pen@epsilon
  pen@params[mask]
}

zapModelHelper <- function(model) {
  if (length(model@submodels)) for (m1 in model@submodels) {
    model <- mxModel(model, zapModelHelper(m1))
  }
  model@penalties <- list()
  model
}

#' mxPenaltyZap
#'
#' Fix any free parameters within \code{epsilon} of zero to
#' zero. These parameters are no longer estimated. Remove all
#' \link[=MxPenalty-class]{MxPenalty} objects from the model.  This is envisioned to be
#' used after using \link{mxPenaltySearch} to locate the best penalty
#' hyperparameters and apply penalties to model estimation. While
#' penalties can simplify a model, they also bias parameters toward
#' zero. By re-estimating the model after using \code{mxPenaltyZap},
#' parameters that remain free are likely to exhibit less bias.
#'
#' @param model an \link{MxModel}
#' @param silent whether to suppress diagnostic output

# alternate names mxPenaltyShrinkToZero mxPenaltyFixSmall
mxPenaltyZap <- function(model, silent=FALSE) {
  penalties <- collectComponents(model, NULL, "penalties", qualifyNames)
  toZap <- Reduce(union, mapply(regularizedToZeroParameters,
                                model@penalties, MoreArgs=list(model)))
  if (!silent && length(toZap)) message(paste("Zapping", omxQuotes(toZap)))
  model <- omxSetParameters(model, labels=toZap, values = 0, free=FALSE)
  toFix <- Reduce(union, mapply(function(pen) pen@hyperparameters, model@penalties))
  if (!silent && length(toFix)) message(paste("Fixing", omxQuotes(toFix)))
  model <- omxSetParameters(model, labels=toFix, free=FALSE)
  if (!silent) message(paste("Tip: Use\n  model = mxRun(model)\nto re-estimate the model",
                             "without any penalty terms."))
  zapModelHelper(model)
}

##' mxPenaltyLASSO
##'
##' Least Absolute Selection and Shrinkage Operator regularization
##'
##' @param lambda strength of the penalty to be applied at starting values (default 0)
##' @param lambda.step step function for lambda step (default .01)
##' @param lambda.max end of lambda range (default .4)
##' @param lambda.min minimum lambda value (default lambda)
##' @template args-dots-barrier
##' @template args-regularize
##'
#' @export
mxPenaltyLASSO <- function(what, name, lambda=0, lambda.step=.01, lambda.max=NA, lambda.min=NA,
                              epsilon=1e-5, scale=1, ..., hyperparams=c('lambda')) {
  prohibitDotdotdot(list(...))
  if(is.na(lambda.min)) lambda.min=lambda
  if(is.na(lambda.max)) lambda.max=lambda.min+40*lambda.step
  if (length(hyperparams) != 1) stop("Provide exactly one hyperparam")
  mxPenalty(what, epsilon, scale, how="lasso", hyperparams=hyperparams, hpranges=list(lambda=seq(lambda, lambda.max, by=lambda.step)), name=name)
}

##' mxPenaltyRidge
##'
##' Ridge regression regularization
##'
##' @param lambda strength of the penalty to be applied at start (default 0)
##' @param lambda.step lambda step during penalty search (default 0.01)
##' @param lambda.max when to end the lambda search (default 0.4)
##' @param lambda.min minimum lambda value (default lambda)
##' @template args-regularize
##' @template args-dots-barrier
##'
#' @export
mxPenaltyRidge <- function(what, name, lambda=0, lambda.step=.01, lambda.max=.4, lambda.min=NA,
                              epsilon=1e-5, scale=1, ..., hyperparams=c('lambda')) {
  prohibitDotdotdot(list(...))
  if(is.na(lambda.min)) lambda.min=lambda
  if (length(hyperparams) != 1) stop("Provide exactly one hyperparam")
  mxPenalty(what, epsilon, scale, how="ridge", hyperparams=hyperparams,
               hpranges = list(lambda=seq(lambda, lambda.max, lambda.step)), name=name)
}

##' mxPenaltyElasticNet
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
##' @template args-dots-barrier
##'
##' @details Applies elastic net regularization.  Elastic net is a weighted combination of ridge and LASSO penalties.
##'
#' @export
mxPenaltyElasticNet <- function(what, name,
                                   alpha=0,  alpha.step=.1,  alpha.max=1,
                                   lambda=0, lambda.step=.1, lambda.max=.4,
                                   alpha.min=NA, lambda.min=NA,
                                   epsilon=1e-5, scale=1, ..., hyperparams=c('alpha', 'lambda')) {
  prohibitDotdotdot(list(...))
  if(is.na(lambda.min)) lambda.min=lambda
  if(is.na(alpha.min)) alpha.min=alpha
  if (length(hyperparams) != 2) stop("Provide exactly two hyperparameters")
  mxPenalty(what, epsilon, scale, how="elasticNet", hyperparams=hyperparams,
               hpranges = list(alpha=seq(alpha, alpha.max, alpha.step), lambda=seq(lambda, lambda.max, by=lambda.step)), name=name)
}

evaluatePenalty <- function(penalty, model, labelsData, env, show, compute) {
  if (!compute) return(penalty@result)

  # MxPenalty is an opaque object like MxFitFunction. show=TRUE should just
  # show the name, not the formula.

  value <- labelToValue(penalty@params, labelsData, model)
  svalue <- abs(value) / penalty@scale
  width <- penalty@epsilon * penalty@smoothProportion
  inactive <- penalty@epsilon - width
  strength <- (svalue - inactive) / width
  strength[strength > 1] <- 1
  strength[strength < 0] <- 0
  lasso <- sum(strength * svalue)
  ridge <- sum(strength * svalue * svalue)

  if (penalty@type == "lasso") {
    lambda <- labelToValue(penalty@hyperparameters[1], labelsData, model)
    result <- lambda * lasso
  } else if (penalty@type == "ridge") {
    lambda <- labelToValue(penalty@hyperparameters[1], labelsData, model)
    result <- lambda * ridge
  } else if (penalty@type == "elasticNet") {
    alpha <- labelToValue(penalty@hyperparameters[1], labelsData, model)
    lambda <- labelToValue(penalty@hyperparameters[2], labelsData, model)
    result <- lambda * ((1-alpha) * ridge + alpha * lasso)
  } else stop(paste("penalty@type", omxQuotes(penalty@type), "is not implemented"))
  result
}

generatePenaltyList <- function(flatModel, modelname, parameters, labelsData) {
	got <- lapply(flatModel@penalties, function(p1) {
    hpnames <- p1@hyperparameters
    lx <- match(p1@hyperparameters, rownames(labelsData))
    if (any(is.na(lx))) {
      stop(paste("Cannot locate penalty parameters",
                 omxQuotes(hpnames[is.na(lx)])))
    }
    coord <- sapply(lx, function(hp1) {
      mat <- match(imxIdentifier(labelsData[hp1, 'model'], labelsData[hp1, 'matrix']),
                   names(flatModel@matrices)) - 1L
      c(mat, labelsData[hp1, 'row'] - 1L, labelsData[hp1, 'col'] - 1L)
    })
    p1@hyperparameters <- c(coord)
    pnames <- p1@params
    p1@params <- match(p1@params, names(parameters)) - 1L
    if (any(is.na(p1@params))) {
      # maybe too strict?
      stop(paste("Cannot locate parameters",
                 omxQuotes(pnames[is.na(p1@params)]), "for penalty search"))
    }
    p1
  })
  if (length(got)) {
    totalParams <- sum(sapply(flatModel@penalties, function(p1) length(p1@params)))
    if (totalParams != length(Reduce(union, lapply(flatModel@penalties, function(p1) p1@params))))
      stop("Attempt to regularize the same parameter more than once")
  }
  got
}
