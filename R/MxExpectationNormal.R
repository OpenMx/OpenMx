#
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

setClass(Class = "BaseExpectationNormal",
         contains = "MxBaseExpectation",
         representation = representation(
           expectedCovariance = "MxOptionalCharOrNumber",
           expectedMean = "MxOptionalCharOrNumber",
           thresholds = "MxCharOrNumber",
           threshnames = "character",
           discrete = "MxCharOrNumber",
           discreteSpec = "MxOptionalMatrix",
           .discreteCheckCount = "logical"))

setMethod("qualifyNames", signature("BaseExpectationNormal"),
          function(.Object, modelname, namespace) {
            for (sl in c(paste0('expected', c('Covariance','Mean')),
                       'thresholds', 'discrete')) {
              if (is.null(slot(.Object, sl))) next
              slot(.Object, sl) <- imxConvertIdentifier(slot(.Object, sl), modelname, namespace)
            }
            .Object
          })

setMethod("genericExpDependencies", signature("BaseExpectationNormal"),
          function(.Object, dependencies) {
            sources <- c(.Object@discrete, .Object@thresholds)
            sources <- sources[!is.na(sources)]
            dependencies <- imxAddDependency(sources, .Object@name, dependencies)
            return(dependencies)
          })

setMethod("genericExpRename", signature("BaseExpectationNormal"),
          function(.Object, oldname, newname) {
            for (sl in c(paste0('expected', c('Covariance','Mean')),
                         'thresholds', 'discrete')) {
              if (is.null(slot(.Object, sl))) next
              slot(.Object,sl) <- renameReference(slot(.Object,sl), oldname, newname)
            }
            .Object
          })

setMethod("genericExpFunConvert", signature("BaseExpectationNormal"),
          function(.Object, flatModel, model, labelsData, dependencies) {
            for (sl in c(paste0('expected', c('Covariance','Mean')))) {
              matName <- slot(.Object,sl)
              if (length(matName) == 0) next
              mat <- flatModel[[matName]]
              if (any(mat@free)) {
                msg <- paste('Free parameters are not allowed in matrix',
                             omxQuotes(matName))
                stop(msg, call.=FALSE)
              }
              if (any(!is.na(mat@labels))) {
                msg <- paste('Labels are not allowed in matrix', omxQuotes(matName))
                stop(msg, call.=FALSE)
              }
            }

            tname <- .Object$thresholds
            dname <- .Object$discrete
            if (!is.na(dname)) {
              dmat <- flatModel[[dname]]
              if (length(colnames(dmat)) == 0) {
                stop(paste(dname, "must have column names"))
              }
              ds <- .Object$discreteSpec
              if (is.null(ds)) stop("discrete parameters but no discreteSpec")
              if (length(colnames(ds)) != length(colnames(dmat))) {
                stop(paste("discreteSpec must have the same number of columns",
                           "as", dname))
              }
              if (any(colnames(ds) != colnames(dmat))) {
                stop(paste("discreteSpec and", dname, "must have the same column names"))
              }
            }
            if (!is.na(tname) && !is.na(dname)) {
              dup <- intersect(colnames(flatModel[[tname]]), colnames(flatModel[[dname]]))
              if (length(dup)) {
                msg <- paste('Manifest', omxQuotes(dup), 'cannot be modelled as',
                             'a threshold and as a discrete variable at the',
                             'same time')
                stop(msg, call.=FALSE)
              }
            }

            .Object
          })

setMethod("genericNameToNumber", signature("BaseExpectationNormal"),
	  function(.Object, flatModel, model) {
		  name <- .Object@name
      for (sl in c(paste0('expected', c('Covariance','Mean')),
                   'thresholds', 'discrete')) {
        slot(.Object,sl) <- imxLocateIndex(flatModel, slot(.Object,sl), name)
      }
		  .Object
	  })

setMethod("genericGetExpected", signature("BaseExpectationNormal"),
	function(.Object, model, what, defvar.row=1, subname=model@name) {
		ret <- list()
		if ('thresholds' %in% what) {
			thrname <- .Object@thresholds
      thr <- matrix( , 0 , 0)
			if(!single.na(thrname)){
				thrname <- .modifyDottedName(subname, thrname, sep=".")
				thr <- mxEvalByName(thrname, model, compute=TRUE, defvar.row=defvar.row)
        tnames <- .Object@threshnames
        if(!single.na(tnames)){
          colnames(thr) <- tnames
        }
			}
      disname <- .Object@discrete
      if (!single.na(disname)) {
        npart <- mxGetExpected(model, c('means','covariance'),
                               defvar.row=defvar.row, subname=subname)
        means <- npart[['means']]
        if (!is.null(means) && ncol(means)==1 && colnames(means)=='one')
          means <- t(means)
        cov <- npart[['covariance']]
        ds <- .Object@discreteSpec
				disname <- .modifyDottedName(subname, disname, sep=".")
				dis <- mxEvalByName(disname, model, compute=TRUE, defvar.row=defvar.row)
        if (length(colnames(ds)) != length(colnames(dis))) {
          stop(paste("discreteSpec must have the same number of columns",
                     "as", disname))
        }
        if (any(colnames(ds) != colnames(dis))) {
          stop(paste("discreteSpec and", disname, "must have the same column names"))
        }

        # discreteSpec:
        #   number of outcomes
        #   1 poisson ; 2 negative binomial prob ; 3 neg binomial mu
        # discrete 3 rows:
        #   zero inflation
        #   poisson: lambda; nbinom: size, prob/mu
        thList <- vector("list", ncol(ds))
        names(thList) <- colnames(ds)
        for (cx in colnames(ds)) {
          if (is.na(ds[1,cx])) {
            stop(paste("discreteSpec[1,",cx,"], the maximum observed count,",
                       "cannot be NA when calcuating the expected distribution"))
          }
          outcome <- 1L:as.integer(ds[1,cx]) - 1L
          zif <- min(max(dis[1,cx],0),1)
          if (ds[2,cx] == 1) {
            pr <- ppois(outcome, dis[2,cx])
          } else if (ds[2,cx] == 2) {
            pr <- pnbinom(outcome, dis[2,cx], prob=dis[3,cx])
          } else if (ds[2,cx] == 3) {
            pr <- pnbinom(outcome, dis[2,cx], mu=dis[3,cx])
          } else { stop(paste("Unknown discrete distribution code", ds[2,cx],
                              "in column", cx)) }
          m1 <- 0
          if (!is.null(means)) m1 <- means[1,cx]
          thList[[cx]] <- p2z(zif + (1-zif) * pr) * sqrt(cov[cx,cx]) + m1
        }
        thBlock <- mxSimplify2Array(thList)
        comb <- matrix(NA,
                       nrow=max(nrow(thr), nrow(thBlock)),
                       ncol=(ncol(thr) + ncol(thBlock)))
        if (length(thr)) {
          comb[1:nrow(thr),1:ncol(thr)] <- thr
          comb[1:nrow(thBlock),ncol(thr) + 1:ncol(thBlock)] <- thBlock
          colnames(comb) <- c(colnames(thr), colnames(thBlock))
          thr <- comb
        } else {
          thr <- thBlock
        }
      }
			ret[['thresholds']] <- thr
      }
    ret
    })

constrainCorData <- function(ex, size, model, flatModel) {
  path <- unlist(strsplit(ex@name, imxSeparatorChar, fixed = TRUE))
  if (path[1] == model@name) {
    submodel <- model
  } else {
    submodel <- model[[ path[1] ]]
  }
  newObj <- model@.newobjects
  # Don't deal with inherited MxData. That's an indication of
  # an advanced user who knows what they are doing.
  if (!is.null(submodel@data) && submodel@data@type == 'cor' && is.null(ex@expectedCovariance)) {
    for (sl in c('expectedCovariance', 'unit_1_by_nObs', 'constraintForCorData')) {
      if (!is.null(submodel[[sl]])) {
        msg <- paste("In model", omxQuotes(submodel@name),
                     "data has type='cor' with", omxQuotes(class(ex)[1]),
                     "but an object named", omxQuotes(sl),
                     "already exists preventing the addition of a constraint",
                     "on the diagonal of the expected covariance")
				stop(msg, call.=FALSE)
      }
    }
    ex@expectedCovariance <- 'expectedCovariance'
    corConstraint <- list(
      ex,
      mxMatrix(name="expectedCovariance", nrow=size, ncol=size),
      mxMatrix(name= "unit_1_by_nObs", type="Unit", ncol=1, nrow=size),
      mxConstraintFromString(name="constraintForCorData", "diag2vec(expectedCovariance) == unit_1_by_nObs"))
    submodel <- mxModel(submodel, corConstraint)
    newObj <- TRUE
  }
  if (path[1] == model@name) {
    model <- submodel
  } else {
    model <- mxModel(model, submodel)
  }
  model@.newobjects <- newObj
  model
}

setClass(Class = "MxExpectationNormal",
	representation = representation(
		covariance = "MxCharOrNumber",
		means = "MxCharOrNumber",
		dims = "character",
		ExpCov = "matrix",
		ExpMean = "matrix",
    numStats = "numeric"),
	contains = "BaseExpectationNormal")

setMethod("initialize", "MxExpectationNormal",
	function(.Object, covariance, means, dims, thresholds, threshnames, discrete,
           discreteSpec, data = as.integer(NA), name = 'expectation') {
		.Object@name <- name
		.Object@covariance <- covariance
		.Object@means <- means
		.Object@data <- data
		.Object@thresholds <- thresholds
		.Object@dims <- dims
		.Object@threshnames <- threshnames
		.Object@discrete <- discrete
		.Object@discreteSpec <- discreteSpec
    .Object@.discreteCheckCount <- TRUE
		.Object@ExpCov <- matrix()
		.Object@ExpMean <- matrix()
		return(.Object)
	}
)

setMethod("qualifyNames", signature("MxExpectationNormal"),
	  function(.Object, modelname, namespace) {
      .Object <- callNextMethod()
		.Object@name <- imxIdentifier(modelname, .Object@name)
		for (s in c("covariance", "means", "data")) {
			if (is.null(slot(.Object, s))) next
			slot(.Object, s) <-
			  imxConvertIdentifier(slot(.Object, s), modelname, namespace)
		}
		return(.Object)
})

setMethod("genericExpDependencies", signature("MxExpectationNormal"),
	function(.Object, dependencies) {
    dependencies <- callNextMethod()
    sources <- c(.Object@covariance, .Object@means)
	sources <- sources[!is.na(sources)]
	dependencies <- imxAddDependency(sources, .Object@name, dependencies)
	return(dependencies)
})

setMethod("genericExpConvertEntities", "MxExpectationNormal",
	function(.Object, flatModel, namespace, labelsData) {
		flatModel <- updateExpectationDimnames(.Object, flatModel, labelsData)
		flatModel <- updateThresholdDimnames(.Object, flatModel, labelsData)
		return(flatModel)
	}
)

setMethod("genericExpRename", signature("MxExpectationNormal"),
	function(.Object, oldname, newname) {
    .Object <- callNextMethod()
		for (s in c("covariance", "means", "data")) {
			if (is.null(slot(.Object, s))) next
			slot(.Object, s) <- renameReference(slot(.Object, s), oldname, newname)
    }
		return(.Object)
})

.modifyDottedName <- function(sub, obj, sep='.'){
	if(single.na(obj)) {return(obj)}
	if(length(obj) > 1) {stop('Please pass me one object name at a time.')}
	if(length(strsplit(obj, split=sep, fixed=TRUE)[[1]]) > 1 ){
		return(obj)
	} else {
		return(paste(sub, obj, sep=sep))
	}
}

setMethod("genericGetExpected", signature("MxExpectationNormal"),
	function(.Object, model, what, defvar.row=1, subname=model@name) {
		ret <- callNextMethod()
		if (any(c('covariance','covariances') %in% what)) {
			covname <- .modifyDottedName(subname, .Object@covariance)
			cov <- mxEvalByName(covname, model, compute=TRUE, defvar.row=defvar.row, .extraBack=7L)
			dnames <- .Object$dims
			if(!single.na(dnames)){
				colnames(cov) <- dnames
				rownames(cov) <- dnames
			}
			ret[['covariance']] <- cov
		}
		if (any(c('means', 'mean') %in% what)) {
			meanname <- .Object@means
			if(!single.na(meanname)){
				meanname <- .modifyDottedName(subname, meanname, sep=".")
				mean <- mxEvalByName(meanname, model, compute=TRUE, defvar.row=defvar.row)
				dnames <- .Object$dims
				if(!single.na(dnames)){
					colnames(mean) <- dnames
				}
			} else {mean <- matrix( , 0, 0)}
			ret[['means']] <- mean
		}
		ret
	})

getThresholdMask <- function(model, cols, subname) {
	if (subname != model$name) {
		d1 <- model[[subname]]@data
	} else {
		d1 <- model@data
	}
	if (is.null(d1)) stop(paste("Cannot find observed thresholds, model",
				    omxQuotes(subname), "has no data"))
	if (d1@type == 'raw') {
		lev <- sapply(d1$observed[,cols], function(x) length(levels(x)))
		mask <- matrix(FALSE, max(lev)-1, length(cols))
		colnames(mask) <- cols
		for (lx in 1:length(cols)) mask[1:(lev[lx]-1), lx] <- TRUE
		mask
	} else {
		th <- d1@thresholds
		if (single.na(th)) stop(paste("Observed data in model", omxQuotes(model$name), "has no thresholds"))
		!is.na(th)
	}
}

setMethod("genericGetExpectedVector", signature("BaseExpectationNormal"),
	function(.Object, model, defvar.row=1, subname=model@name) {
		ret <- genericGetExpected(.Object, model, c('covariance', 'means', 'thresholds'), defvar.row, subname)
		cov <- ret[['covariance']]
		mns <- ret[['means']]
		if (is.null(mns)) stop("mns is null")
		nv <- nrow(cov)
		covNames <- paste0('cov', vech(outer(1:nv, 1:nv, FUN=paste, sep='_')))
		mnsNames <- paste0('mean', 1:nv)
		thr <- ret[['thresholds']]
		if (prod(dim(thr)) == 0) {
			v <- c(vech(cov), mns[!is.na(mns)])
			names(v) <- c(covNames, mnsNames[!is.na(mns)])
		} else {
			thrNames <- outer(paste0('thr', 1:nrow(thr)), 1:ncol(thr), paste, sep='_')
			dth <- getThresholdMask(model, colnames(thr), subname)
			v <- c(vech(cov), mns[!is.na(mns)], thr[dth])
			names(v) <- c(covNames, mnsNames[!is.na(mns)], thrNames[dth])
		}
		return(v)
})

setMethod("genericGetExpectedStandVector", signature("BaseExpectationNormal"),
	function(.Object, model, defvar.row=1, subname=model@name) {
		ret <- genericGetExpected(.Object, model, c('covariance', 'means', 'thresholds'), defvar.row, subname)
		cov <- ret[['covariance']]
		mns <- ret[['means']]
		if (is.null(mns)) stop("mns is null")
		nv <- nrow(cov)
		thr <- ret[['thresholds']]
		if (prod(dim(thr)) != 0) {
			thrNames <- outer(paste0('thr', 1:nrow(thr)), 1:ncol(thr), paste, sep='_')
			dth <- getThresholdMask(model, colnames(thr), subname)
		}
		v <- .standardizeCovMeansThresholds(cov, mns, thr, dth, vector=TRUE)
		return(v)
})

.standardizeCovMeansThresholds <- function(cov, means, thresholds, dth, vector=FALSE){
  mnames <- colnames(cov)
  if (length(mnames) == 0) {
    stop("I give up. Have no idea how to standardize this expectation.\nYour covariance matrix must have dimnames.")
  }
  if (prod(dim(thresholds)) == 0) {
    ordInd <- c()
  } else {
    ordInd <- match(colnames(thresholds), mnames)
    thresholds <- matrix( (c(thresholds) - rep(means[ordInd], each=nrow(thresholds)) ) / rep(sqrt(diag(cov)[ordInd]), each=nrow(thresholds)), nrow=nrow(thresholds), ncol=ncol(thresholds) )
    means[ordInd] <- means[ordInd] - means[ordInd]
    cov <- .ordinalCov2Cor(cov, ordInd)
  }
	if(!vector){
		return(list(cov=cov, means=means, thresholds=thresholds))
	} else {
		v <- c()
		vn <- c()
		if (length(means)) {
		  for (vx in 1:length(mnames)) {
		    tcol <- which(vx == ordInd)
		    if (length(tcol) == 0) {
		      v <- c(v, means[vx])
		      vn <- c(vn, mnames[vx])
		    } else {
		      tcount <- sum(dth[,tcol])
		      v <- c(v, thresholds[1:tcount,tcol])
		      vn <- c(vn, paste0(mnames[vx], 't', 1:tcount))
		    }
		  }
		}
		for (vx in 1:length(mnames)) {
			if (any(vx == ordInd)) next
			v <- c(v, cov[vx,vx])
			vn <- c(vn, paste0('var_', mnames[vx]))
		}
		v <- c(v, vechs(cov))
		nv <- length(mnames)
		vn <- c(vn, paste0('poly_', vechs(outer(mnames[1:nv], mnames[1:nv], FUN=paste, sep='_'))))
		names(v) <- vn
		return(v)
	}
}

.ordinalCov2Cor <- function(cov, ordInd){
	dim <- ncol(cov)
	egOutCov <- matrix(0, nrow=dim, ncol=dim)
	stddev <- sqrt(diag(cov))
	if(is.logical(ordInd)){notOrdInd <- !ordInd} else {notOrdInd <- -ordInd}
	stddev[notOrdInd] <- 1
	for(i in 1:dim) {
		for(j in 1:i) {
			egOutCov[i,j] = cov[i, j] / (stddev[i] * stddev[j]);
			egOutCov[j,i] = egOutCov[i,j]
		}
	}
	diag(egOutCov)[ordInd] <- 1
	return(egOutCov)
}

imxGetExpectationComponent <- function(model, component, defvar.row=1, subname=model$name)
{
  warnModelCreatedByOldVersion(model)
  if (!is(model[[subname]], "MxModel")) {
    stop(paste("Submodel", subname, "in model", model$name, "is not found"))
  }
	if(is.null(model[[subname]]$expectation) && (class(model[[subname]]$fitfunction) %in% "MxFitFunctionMultigroup") ){
		submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
		got <- list()
		for(amod in submNames){
			got[[amod]] <- imxGetExpectationComponent(model, component, defvar.row=1, subname=amod)
		}
		if(component=='vector' || tolower(component)=='standvector'){got <- unlist(got)}
		got
	} else if (length(component) == 1 && component == 'vector') {
		genericGetExpectedVector(model[[subname]]$expectation, model, defvar.row, subname)
	} else if (length(component) == 1 && tolower(component) == 'standvector') {
		genericGetExpectedStandVector(model[[subname]]$expectation, model, defvar.row, subname)
	} else {
		got <- genericGetExpected(model[[subname]]$expectation, model, component, defvar.row, subname)
		if (length(got) == 1) {
			got[[1]]
		} else {
			got
		}
	}
}

mxGetExpected <- imxGetExpectationComponent

sse <- function(x){sum(x^2)}

#' Estimate the Jacobian of manifest model with respect to parameters
#'
#' The manifest model excludes any latent variables or processes. For
#' RAM and LISREL models, the manifest model contains only the
#' manifest variables with free means, covariance, and thresholds.
#'
#' @details
#' The Jacobian is estimated by the central finite difference.
#'
#' If the \code{standardize} argument is TRUE, then the Jacobian is for the standardized model.
#' For Normal expectations the standardized manifest model has the covariances returned as correlations, the variances returned as ones, the means returned as zeros, and the thresholds are returned as z-scores.
#' For the thresholds the z-scores are computed by using the model-implied means and variances.
#'
#' @param model an mxModel
#' @param defvar.row which row to use for definition variables
#' @param standardize logical, whether or not to standardize the parameters
#' @return a matrix with manifests in the rows and original parameters in the columns
#' @seealso \link{mxGetExpected}
omxManifestModelByParameterJacobian <- function(model, defvar.row=1, standardize=FALSE) {
	theParams <- omxGetParameters(model)
	if (standardize) {
		ex <- 'expectation'
		if (is.null(model$expectation) && is(model$fitfunction, "MxFitFunctionMultigroup")) {
			submNames <- sapply(strsplit(model$fitfunction$groups, ".", fixed=TRUE), "[", 1)
			ex <- paste0(submNames, ".expectation")
		}
		tmpModel <- mxModel(model, mxComputeJacobian(defvar.row=defvar.row, of=ex))
		tmpModel <- mxRun(tmpModel, silent=TRUE)
		jac <- tmpModel$compute$output$jacobian
		dimnames(jac) <- list(names(mxGetExpected(model, 'standVector')), names(theParams))
	} else {
		jac <- numDeriv::jacobian(func=.mat2param, x=theParams, method.args=list(r=2), model=model, defvar.row=defvar.row, standardize=standardize)
		dimnames(jac) <- list(names(mxGetExpected(model, 'vector')), names(theParams))
	}

	return(jac)
}

mxCheckIdentification <- function(model, details=TRUE){
  warnModelCreatedByOldVersion(model)
	notAllowedFits <- c("MxFitFunctionAlgebra", "MxFitFunctionRow", "MxFitFunctionR")
	if( class(model$fitfunction) %in% notAllowedFits ){
		msg <- paste("Identification check is not possible for models with", omxQuotes(notAllowedFits), 'fit functions.\n', "If you have a multigroup model, use mxFitFunctionMultigroup.")
		stop(msg, call.=FALSE)
	}
	if(imxHasDefinitionVariable(model)){
		stop("Beep beep ribby ribby.  I found definition variables in your model.\nI might not give you the identification answer you're looking for in this case.  See ?mxCheckIdentification.")
	}
	eps <- 1e-17
	theParams <- omxGetParameters(model)
	jac <- omxManifestModelByParameterJacobian(model)
	if(imxHasConstraint(model)){
		tmpModel <- mxModel(model, mxComputeSequence(list(mxComputeNumericDeriv(hessian=FALSE), mxComputeReportDeriv())))
		tmpModel <- mxRun(tmpModel, silent=TRUE)
		cjac <- tmpModel$output$constraintJacobian
		# drop model name from constraint name
		rownames(cjac) <- sapply(strsplit(rownames(cjac), fixed=TRUE, split=imxSeparatorChar), '[', 2)
	} else {
		cjac <- matrix(, nrow=0, ncol=length(theParams))
		colnames(cjac) <- names(theParams)
	}
	# Concatenate via rows the model and constraint Jacobians
	jac <- rbind(jac, cjac)
	# Check that rank of jac == length(theParams)
	rank <- qr(jac)$rank
	if(rank == length(theParams)){
		message("Model is locally identified")
		stat <- TRUE
	} else {
		message("Model is not locally identified")
		stat <- FALSE
	}
	if(details == TRUE){
		jacOC <- Null(t(jac)) # Orthogonal complement of t(jac), i.e. the basis for the null space of the column space of jac
		nidp <- names(theParams)[apply(jacOC, 1, sse) > eps] # non-identified free params have non-zero rows
		if(length(nidp) == 0) {
			nidp <- "None"
		}
	} else {
		nidp <- "Not Requested"
	}
	return(list(status=stat, jacobian=jac, non_identified_parameters=nidp))
}

.mat2param <- function(x, model, defvar.row=1, standardize=FALSE){
  paramNames <- names(omxGetParameters(model))
  model <- omxSetParameters(model, values=x, labels=paramNames, free=TRUE)
  if(!standardize){
    got <- mxGetExpected(model, 'vector', defvar.row)
  } else {
    got <- mxGetExpected(model, 'standVector', defvar.row)
  }
  return(got)
}

setGeneric("genericGenerateData",
	function(.Object, model, nrows, subname, empirical, returnModel, use.miss,
		   .backend, nrowsProportion) {
	return(standardGeneric("genericGenerateData"))
})

setMethod("genericGenerateData", signature("MxExpectationNormal"),
	  function(.Object, model, nrows, subname, empirical, returnModel, use.miss,
		   .backend, nrowsProportion) {
	    return(generateNormalData(model, nrows, subname, empirical, returnModel,
				      use.miss, .backend, nrowsProportion))
	  })

.rmvnorm <- function(nrow, mean, sigma, empirical) {
  if (!empirical) {
    mvtnorm::rmvnorm(nrow, mean, sigma)
  } else {
    MASS::mvrnorm(nrow, mu=mean, Sigma=sigma, empirical=TRUE)
  }
}

calcNumRows <- function(nrows, nrowsProportion, origRows, subname) {
  if (is.null(nrows)) {
    if (is.null(origRows)) stop("You must specify nrows")
    if (!is.null(nrowsProportion)) {
      nrows <- round(origRows * nrowsProportion)
      if (nrows < 1) {
	nrows <- 1
	warning(paste("In model", omxQuotes(subname),
		      "nrowsProportion is too small. Rounding up to 1 row"))
      }
    } else {
      nrows <- origRows
    }
  }
  nrows
}

generateNormalData <- function(model, nrows, subname, empirical, returnModel, use.miss,
			       .backend, nrowsProportion) {
  origData <- findDataForSubmodel(model, subname)
  origRows <- if (!is.null(origData)) { nrowMxData(origData) } else { NULL }
  nrows <- calcNumRows(nrows, nrowsProportion, origRows, subname)

	# Check for definition variables
	if(imxHasDefinitionVariable(model[[subname]])){
		if (origData$type != 'raw') {
			stop(paste("Definition variable(s) found, but original data is type",
				omxQuotes(origData$type)))
		}
		origData <- origData$observed
		if(nrows != nrow(origData)){
			stop("Definition variable(s) found, but the number of rows in the data do not match the number of rows requested for data generation.")
		}
		# Generate data row by row
		theCov <- imxGetExpectationComponent(model, "covariance", subname=subname)
		data <- matrix(NA, nrow=nrows, ncol=ncol(theCov))
		colnames(data) <- colnames(theCov)
		data <- as.data.frame(data)
		for(i in 1:nrows){
			theMeans <- imxGetExpectationComponent(model, "means", defvar.row=i, subname=subname)
			theCov <- imxGetExpectationComponent(model, "covariance", defvar.row=i, subname=subname)
			theThresh <- imxGetExpectationComponent(model, "thresholds", defvar.row=i, subname=subname)
			data[i,] <- .rmvnorm(1, theMeans, theCov, empirical)
		}
		data <- ordinalizeDataHelper(data, theThresh, origData)
		if (!is.null(origData)) {
			for (dcol in setdiff(colnames(origData), colnames(data))) {
				data[[dcol]] <- origData[[dcol]]
			}
		}
	} else if (!is.null(origData) && origData$type == 'cov' && returnModel) {
		theMeans <- imxGetExpectationComponent(model, "means", subname=subname)
		theCov <- imxGetExpectationComponent(model, "covariance", subname=subname)
		if (length(theMeans) == 0) theMeans <- NA
		if (empirical) {
		  got <- mxModel(model[[subname]], mxData(theCov, "cov", means=theMeans, numObs=nrows))
		} else {
		  newCov <- rWishart(1, nrows, theCov)[,,1] / nrows
		  dimnames(newCov) <- dimnames(theCov)
		  got <- mxModel(model[[subname]],
				 mxData(newCov, "cov", means=theMeans, numObs=nrows))
		}
		  return(got)
	} else{
		theMeans <- imxGetExpectationComponent(model, "means", subname=subname)
		theCov <- imxGetExpectationComponent(model, "covariance", subname=subname)
		theThresh <- imxGetExpectationComponent(model, "thresholds", subname=subname)
		if (length(theMeans) == 0) {
			theMeans <- rep(0, nrow(theCov))
		}
		data <- .rmvnorm(nrows, theMeans, theCov, empirical)
		colnames(data) <- colnames(theCov)
		data <- as.data.frame(data)
		data <- ordinalizeDataHelper(data, theThresh, origData)
	}
  if (use.miss && !is.null(origData) && origData$type == 'raw' &&
        all(colnames(data) %in% colnames(origData$observed))) {
	  del <- is.na(origData$observed[,colnames(data),drop=FALSE])
	  if (nrows != origRows) {
	    del    <- del[sample.int(origRows, nrows, replace=TRUE),,drop=FALSE]
	  }
	  data[del] <- NA
	}
	if (returnModel) {
	  mxModel(model[[subname]], mxData(as.data.frame(data), "raw"))
	} else {
	  as.data.frame(data)
	}
}

ordinalizeDataHelper <- function(data, thresh, origData=NULL) {
	if (!is.null(origData) && is(origData, "MxData")) origData <- origData$observed
	if( prod(dim(thresh)) != 0){
		ordvars <- colnames(thresh)
    if (length(ordvars) != ncol(thresh)) stop("thresholds missing column names")
		for(avar in ordvars){
			delthr <- thresh[,avar]
			usethr <- 1:sum(!is.na(delthr))  # assumes NA indicates unused threshold
			if (!is.null(origData) && !is.null(origData[[avar]])) {
				usethr <- 1:(length(levels(origData[[avar]])) - 1L)
			}
			delthr <- delthr[usethr]
			levthr <- 1L:(length(usethr)+1L)
			if (!is.null(origData)  && !is.null(origData[[avar]])) {
				levthr <- levels(origData[[avar]])
			}
			delvar <- cut(as.vector(data[,avar]), c(-Inf, delthr, Inf), labels=levthr)
			data[,avar] <- mxFactor(delvar, levels=levthr)
		}
	}
	return(data)
}

generateRelationalData <- function(model, returnModel, .backend, subname, empirical) {
	model <- model[[subname]]
	# TODO add outside data reference to relational data
	if (.backend) {
	  if (empirical) stop("mxGenerateData(..., empirical=TRUE) is not implemented for multilevel models")

		plan <- mxComputeGenerateData()
		model$expectation$.maxDebugGroups <- 0L
		model$expectation$.optimizeMean <- 0L
		modelE <- mxModel(model, plan)
		modelE <- mxRun(modelE, silent=TRUE)
		simData <- modelE$compute$output

		datalist <- modelE@runstate$datalist
		for (dName in names(datalist)) {
			if (is.null(simData[[dName]])) {
				simData[[dName]] <- datalist[[dName]]$observed
			} else {
				orig <- datalist[[dName]]$observed
				toCopy <- setdiff(colnames(orig), colnames(simData[[dName]]))
				for (col in toCopy) {
					simData[[dName]][[col]] <- orig[[col]]
				}
			}
		}

		names(simData) <- substr(names(simData), 1, nchar(names(simData))-5) #strip .data

		if (!returnModel) {
			return(simData)
		} else {
			for (modelName in names(simData)) {
				if (modelName == model$name) {
					model@data@observed <- simData[[modelName]]
				} else {
					model[[modelName]]@data@observed <- simData[[modelName]]
				}
			}
			return(model)
		}
	}

	# nocov start
	plan <- mxComputeSequence(list(
	    mxComputeOnce('expectation', 'distribution', 'flat'),
	    mxComputeReportExpectation()
	))

	modelE <- mxModel(model, plan)
	modelE$expectation$.rampartCycleLimit <- 0L
	modelE <- mxRun(modelE, silent=TRUE)
	dataEnv <- new.env()
	for (dName in names(modelE@runstate$datalist)) {
		modelName <- substr(dName, 1, nchar(dName)-5)  # remove .data
		assign(modelName, modelE@runstate$datalist[[ dName ]]@observed, envir=dataEnv)
	}
	ed <- modelE$expectation$debug
	layout <- ed$layout
	fmt <- paste0('g%0', ceiling(log10(ed$numGroups)), 'd')
	for (gx in 1:ed$numGroups) {
		groupName <- sprintf(fmt, gx)
		clumpSize <- ed[[groupName]]$clumpSize
		numCopies <- nrow(layout[layout$group == gx,]) %/% clumpSize
		cxLength <- length(ed[[groupName]]$mean) %/% numCopies
		groupTodo <- ed$layout[ed[[groupName]]$layout[,'aIndex'],]
		for (cx in 1:numCopies) {
			todo <- groupTodo[seq(1+(cx-1)*clumpSize, cx*clumpSize),]
			repl1 <- .rmvnorm(1, ed[[groupName]]$mean[seq(1+(cx-1)*cxLength, cx*cxLength)],
					  sigma=as.matrix(ed[[groupName]]$covariance), empirical)
			dx <- 1
			for (tx in 1:nrow(todo)) {
				modelName <- as.character(todo[tx,'model'])
				if (modelName == modelE$name) {
					submodel <- modelE
				} else {
					submodel <- modelE[[modelName]]
				}
				row <- todo[tx,'row']
				manifests <- rownames(submodel$F)
				beforeData <- dataEnv[[modelName]][row, manifests]
				notMissing <- !is.na(beforeData)
				if (sum(notMissing) > 0) {
					afterData <- repl1[seq(dx, dx+sum(notMissing) - 1)]
					dataEnv[[modelName]][row, manifests[notMissing] ] <- afterData
					dx <- dx + sum(notMissing)
				}
			}
		}
	}
	if (!returnModel) {
		ret <- list()
		for (n in names(dataEnv)) {
			ret[[n]] <- dataEnv[[n]]
		}
		ret
	} else {
		for (modelName in names(dataEnv)) {
			if (modelName == model$name) {
				model@data@observed <- dataEnv[[modelName]]
			} else {
				model[[modelName]]@data@observed <- dataEnv[[modelName]]
			}
		}
		model
	}
	# nocov end
}

simulate.MxModel <- function(object, nsim = 1, seed = NULL, ...) {
	if (!is.null(seed)) {
		set.seed(seed)
	}
	mxGenerateData(object, nsim)
}

extractData <- function(model) {
	datasets <- list()
	if (!is.null(model@data)) datasets <- c(datasets, list(model@data))
	if (length(model@submodels)) {
		datasets <- c(datasets, unlist(lapply(model@submodels, extractData), recursive=FALSE))
	}
	return(datasets)
}

extractObservedData <- function(model) {
	lapply(extractData(model), function(mxd) {
		if (mxd@type != 'raw') stop(paste("Cannot extract observed data when type=", mxd@type))
		mxd@observed
	})
}

mxGenerateData <- function(model, nrows=NULL, returnModel=FALSE, use.miss = TRUE,
			   ..., .backend=TRUE, subname=NULL, empirical=FALSE, nrowsProportion=NULL) {
	prohibitDotdotdot(list(...))
	if (!is.null(nrows) && !is.null(nrowsProportion)) {
	  stop("You cannot specify both nrows and nrowsProportion")
	}
	if (!is.null(nrows) && nrows != round(nrows)) {
	  stop(paste("Cannot generate a non-integral number of rows:", nrows))
	}
	if (is(model, 'data.frame')) {
		fake <- omxAugmentDataWithWLSSummary(mxData(model,'raw'), "ULS", "marginals",
			fullWeight=FALSE, returnModel=TRUE)
		obsStats <- fake$data$observedStats
		fake$S$values <- obsStats$cov
		fake$M$values <- obsStats$means
		if (!is.null(obsStats$thresholds)) fake$thresh$values <- obsStats$thresholds

		if (!is.null(nrowsProportion)) {
		  nrows <- round(nrow(model) * nrowsProportion)
		  if (nrows < 1) {
		    nrows <- 1
		    warning(paste("In model", omxQuotes(model$name),
				  "nrowsProportion is too small. Rounding up to 1 row"))
		  }
		}
		if(is.null(nrows)) nrows <- nrow(model)
		return(mxGenerateData(fake, nrows, returnModel, empirical=empirical))
	}
  warnModelCreatedByOldVersion(model)
	if(is.null(subname)) subname <- model$name
	if (is.null(model[[subname]]$expectation) && is(model[[subname]]$fitfunction, 'MxFitFunctionMultigroup')) {
		todo <- sub(".fitfunction", "", model[[subname]]$fitfunction$groups, fixed=TRUE)
		for (s1 in todo) {
		  model <- mxModel(model, mxGenerateData(model, returnModel=TRUE, nrows=nrows,
							 use.miss=use.miss, .backend=.backend, subname=s1,
							 empirical=empirical, nrowsProportion=nrowsProportion))
		}
		if (!returnModel) {
			return(extractObservedData(model))
		} else {
			return(model)
		}
	}
	genericGenerateData(model[[subname]]$expectation, model, nrows, subname, empirical,
			    returnModel, use.miss, .backend, nrowsProportion)
}

verifyExpectedObservedNames <- function(data, covName, flatModel, modelname, objectiveName) {
	covariance <- flatModel[[covName]]
	if (is(covariance, "MxMatrix") && !identical(dim(covariance), dim(data))) {
		msg <- paste("The dimensions for the expected covariance matrix",
			"and the observed covariance matrix",
			"in the", objectiveName, "expectation function in model",
			omxQuotes(modelname), "are not identical.")
		stop(msg, call. = FALSE)
	}
	if (!identical(dimnames(covariance), dimnames(data))) {
		msg <- paste0("The dimnames for the expected covariance matrix (", omxQuotes(rownames(covariance)),
			") and the observed covariance matrix (", omxQuotes(rownames(data)),
			") in the ", objectiveName, " expectation function in model ",
			omxQuotes(modelname), " are not identical.")
		stop(msg, call. = FALSE)
	}
}

verifyMeans <- function(meansName, mxDataObject, flatModel, modelname) {
	means <- flatModel[[meansName]]
	if (!is.null(means)) {
		if(any(is.na(mxDataObject@means))) {
			msg <- paste("In model", omxQuotes(modelname),
				"the Normal expectation function contains an expected means",
				"vector but the model is missing some data",
				"for the observed means.")
			stop(msg, call. = FALSE)
		}
		meanDimnames <- dimnames(means)
	}
}

setMethod("genericExpFunConvert", "MxExpectationNormal",
	function(.Object, flatModel, model, labelsData, dependencies) {
    .Object <- callNextMethod()
		modelname <- imxReverseIdentifier(model, .Object@name)[[1]]
		name <- .Object@name
		if(is.na(.Object@data)) {
			msg <- paste("The normal expectation function",
				"does not have a dataset associated with it in model",
				omxQuotes(modelname))
			stop(msg, call.=FALSE)
		}

		mxDataObject <- flatModel@datasets[[.Object@data]]
		dataName <- .Object@data
		.Object@data <- imxLocateIndex(flatModel, .Object@data, name)
		threshName <- .Object@thresholds
		covName <- .Object@covariance
		covariance <- flatModel[[covName]]
		.Object@covariance <- imxLocateIndex(flatModel, .Object@covariance, name)
		meansName <- .Object@means
		.Object@means <- imxLocateIndex(flatModel, .Object@means, name)

		if (inherits(mxDataObject, "MxDataDynamic")) return(.Object)

		if (mxDataObject@type != "raw") {
			verifyExpectedObservedNames(mxDataObject@observed, covName, flatModel, modelname, "Normal")
			verifyMeans(meansName, mxDataObject, flatModel, modelname)
		}
		verifyObservedNames(mxDataObject@observed, mxDataObject@means, mxDataObject@type, flatModel, modelname, "Normal")
		checkNumericData(mxDataObject)
		covNames <- colnames(covariance)
		verifyMvnNames(covName, meansName, "expected", flatModel, modelname, class(.Object))
		.Object@dataColumnNames <- covNames
		.Object@dataColumns <- generateDataColumns(flatModel, covNames, dataName)
		verifyThresholds(flatModel, model, labelsData, dataName, covNames, threshName)
		if (single.na(.Object@dims)) {
			.Object@dims <- covNames
		}
		return(.Object)
})

verifyMvnNames <- function(covName, meansName, type, flatModel, modelname, expectationName) {
	if (is.na(meansName)) {
		means <- NA
	} else {
		means <- flatModel[[meansName]]
	}
	covariance <- flatModel[[covName]]
	if (length(covariance)) {
		covDimnames <- dimnames(covariance)
		if (is.null(covDimnames)) {
			msg <- paste("The",type,"covariance matrix associated",
				     "with", expectationName, "in model",
				     omxQuotes(modelname), "does not contain dimnames.")
			stop(msg, call. = FALSE)
		}
		covRows <- covDimnames[[1]]
		covCols <- covDimnames[[2]]
		if (is.null(covRows) || is.null(covCols) ||
		    (length(covRows) != length(covCols)) || !all(covRows == covCols)) {
			msg <- paste("The",type,"covariance matrix associated",
				     "with", expectationName, "in model",
				     omxQuotes(modelname), "does not contain identical",
				     "row and column dimnames.")
			stop(msg, call.=FALSE)
		}
	}
	if (is.null(means) || (!isS4(means) && is.na(means)) || !length(means)) return()
	meanDimnames <- dimnames(means)
	if (is.null(meanDimnames)) {
			msg <- paste("The",type,"means matrix associated",
				"with", expectationName, "in model",
				omxQuotes(modelname), "does not contain dimnames.")
			stop(msg, call.=FALSE)
	}
	meanRows <- meanDimnames[[1]]
	meanCols <- meanDimnames[[2]]
	meanNames <- c()
	if (length(meanRows) > length(meanCols)) { meanNames <- meanRows } else { meanNames <- meanCols }
	if ((length(covCols) != length(meanNames)) || !all(covCols == meanNames)) {
			msg <- paste("The",type,"covariance and",type,
				"means matrices associated",
				"with", expectationName, "in model",
				omxQuotes(modelname), "do not contain identical",
				"dimnames.")
			stop(msg, call.=FALSE)
	}
}

verifyObservedNames <- function(data, means, type, flatModel, modelname, expectationName) {
	dataNames <- dimnames(data)
	if(is.null(dataNames)) {
		msg <- paste("The observed data associated with the",
			expectationName, "expectation function in model",
			omxQuotes(modelname), "does not contain dimnames.")
		stop(msg, call. = FALSE)
	}
	if (type == "cov" || type == "cor") {
		if (length(dataNames) < 2 ||
			is.null(dataNames[[1]]) || is.null(dataNames[[2]]) ||
			!identical(dataNames[[1]], dataNames[[2]])) {
				msg <- paste("The dataset associated with the", expectationName,
					"expectation function in model", omxQuotes(modelname),
    	            "does not contain identical row and column non-NULL dimnames.")
			stop(msg, call. = FALSE)
		}
		if (!single.na(means) && is.null(dimnames(means))) {
			msg <- paste("In model", omxQuotes(modelname),
				", the observed means vector does not contain column names.",
				"Use the names() function to assign names to the means vector.")
			stop(msg, call. = FALSE)
		}
		if (!single.na(means) && !identical(dataNames[[1]], dimnames(means)[[2]])) {
			msg <- paste("The observed covariance or correlation matrix associated with the", expectationName,
				"expectation function in model", omxQuotes(modelname),
				"does not contain identical dimnames to the observed means vector.")
			stop(msg, call. = FALSE)
		}
	} else if ((type == "raw") && (length(dataNames) < 2 || is.null(dataNames[[2]]))) {
		msg <- paste("The dataset associated with the", expectationName,
				"expectation function in model", omxQuotes(modelname),
				"does not contain column names (use dimnames).")
		stop(msg, call. = FALSE)
	}
}

generateDataColumns <- function(flatModel, covNames, dataName) {
	retval <- c()
	if (length(covNames) == 0) return(retval)
	dataColumnNames <- colnames(flatModel@datasets[[dataName]]@observed)
	retval <- match(covNames, dataColumnNames)
	if (any(is.na(retval))) {
		msg <- paste("The column name(s)", omxQuotes(covNames[is.na(retval)]),
			     "in the expected covariance matrix",
			     "of the expectation function",
			     "cannot be found in the column names of",
			     omxQuotes(dataName))
		stop(msg, call. = FALSE)
	}
	return(retval - 1L)
}


updateThresholdDimnames <- function(flatExpectation, flatModel, labelsData) {
	threshName <- flatExpectation@thresholds
	if (is.na(threshName)) {
		return(flatModel)
	}
	thresholds <- flatModel[[threshName]]
	if (is.null(thresholds)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown thresholds name",
			omxQuotes(simplifyName(threshName, modelname)),
			"detected in the expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatExpectation@threshnames
	if (!is.null(colnames(thresholds)) && !single.na(dims) &&
		!identical(colnames(thresholds), dims)) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The thresholds matrix associated",
		"with the expectation function in model",
		omxQuotes(modelname), "contains column names and",
		"the expectation function has specified non-identical threshnames.")
		stop(msg, call.=FALSE)
	}
	if (is.null(colnames(thresholds)) && !single.na(dims)) {
		if (!flatModel@unsafe) {
			tuple <- evaluateMxObject(threshName, flatModel, labelsData, new.env(parent = emptyenv()))
			threshMatrix <- tuple[[1]]
			if (ncol(threshMatrix) != length(dims)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The thresholds matrix associated",
					"with the expectation function in model",
					omxQuotes(modelname), "is not of the same length as the 'threshnames'",
					"argument provided by the expectation function. The 'threshnames' argument is",
					"of length", length(dims), "and the expected covariance matrix",
					"has", ncol(threshMatrix), "columns.")
				stop(msg, call.=FALSE)
			}
		}
		dimnames(flatModel[[threshName]]) <- list(NULL, dims)
	}
	return(flatModel)
}

updateExpectationDimnames <- function(flatExpectation, flatModel,
		labelsData) {
	unsafe <- flatModel@unsafe
	covName <- flatExpectation@covariance
	meansName <- flatExpectation@means
	if (is.na(meansName)) {
		means <- NA
	} else {
		means <- flatModel[[meansName]]
	}
	covariance <- flatModel[[covName]]
	if (is.null(covariance)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown expected covariance name",
			omxQuotes(simplifyName(covName, modelname)),
			"detected in the expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	if (is.null(means)) {
		modelname <- getModelName(flatExpectation)
		stop(paste("Unknown expected means name",
			omxQuotes(simplifyName(meansName, modelname)),
			"detected in the expectation function",
			"of model", omxQuotes(modelname)), call. = FALSE)
	}
	dims <- flatExpectation@dims
	if (!is.null(dimnames(covariance)) && !single.na(dims) &&
		!identical(dimnames(covariance), list(dims, dims))) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The expected covariance matrix associated",
			"with the expectation function in model",
			omxQuotes(modelname), "contains dimnames: ",
            paste(toString(dimnames(covariance)), ".", sep = ""),
			"The expectation function has specified dimnames:",
			paste(toString(dims), ".", sep =""))
		stop(msg, call.=FALSE)
	}
	if (is.null(dimnames(covariance)) && !single.na(dims)) {
		if (!unsafe) {
			tuple <- evaluateMxObject(covName, flatModel, labelsData, new.env(parent = emptyenv()))
			covMatrix <- tuple[[1]]
			if (nrow(covMatrix) != ncol(covMatrix)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected covariance matrix associated",
					"with the expectation function in model",
					omxQuotes(modelname), "is not a square matrix.")
				stop(msg, call.=FALSE)
			}
			if (nrow(covMatrix) != length(dims)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected covariance matrix associated",
					"with the expectation function in model",
					omxQuotes(modelname), "is not of the same length as the 'dimnames'",
					"argument provided by the expectation function. The 'dimnames' argument is",
					"of length", length(dims), "and the expected covariance matrix",
					"has", nrow(covMatrix), "rows and columns.")
				stop(msg, call.=FALSE)
			}
		}
		dimnames(flatModel[[covName]]) <- list(dims, dims)
	}

	if (!isS4(means) && is.na(means)) {
		return(flatModel)
	}

	if (!is.null(dimnames(means)) && !single.na(dims) &&
		!identical(dimnames(means), list(NULL, dims))) {
		modelname <- getModelName(flatExpectation)
		msg <- paste("The expected means matrix associated",
			"with the expectation function in model",
			omxQuotes(modelname), "contains dimnames: ",
            paste(toString(dimnames(means)), ".", sep = ""),
			"The expectation function has specified dimnames:",
			paste(toString(dims), ".", sep =""))
		stop(msg, call.=FALSE)
	}
	if (is.null(dimnames(means)) && !single.na(dims)) {
		if (!unsafe) {
			tuple <- evaluateMxObject(meansName, flatModel, labelsData, new.env(parent = emptyenv()))
			meansMatrix <- tuple[[1]]
			if (nrow(meansMatrix) != 1) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected means vector associated",
					"with the expectation function in model",
					omxQuotes(modelname), "is not a 1 x n matrix.",
					"It has dimensions", nrow(meansMatrix), "x",
					paste(ncol(meansMatrix), '.', sep=''))
				stop(msg, call.=FALSE)
			}
			if (ncol(meansMatrix) != length(dims)) {
				modelname <- getModelName(flatExpectation)
				msg <- paste("The expected means vector associated",
					"with the expectation function in model",
					omxQuotes(modelname), "is not of the same length as the 'dimnames'",
					"argument provided by the expectation function. The 'dimnames' argument is",
					"of length", length(dims), "and the expected means vector",
					"has", ncol(meansMatrix), "columns.")
				stop(msg, call.=FALSE)
			}
		}
		dimnames(flatModel[[meansName]]) <- list(NULL, dims)
	}
	return(flatModel)
}

checkThreshnames <- function(threshnames) {
	if (single.na(threshnames)) threshnames <- as.character(NA)
	if (!is.vector(threshnames) || typeof(threshnames) != 'character') {
		stop("'threshnames' argument is not a character vector")
	}
	if (length(threshnames) == 0) {
		stop("'threshnames' argument cannot be an empty vector")
	}
	if (length(threshnames) > 1 && any(is.na(threshnames))) {
		stop("NA values are not allowed for 'threshnames' vector")
	}
	tt <- table(threshnames)
	if (any(tt > 1)) {
		stop(paste("'threshnames' argument contains", omxQuotes(names(tt)[tt > 1]),
			   "more than once. \nIf you are having problems with Doppelgangers",
			   "perhaps you should check the basement for pods :)"))
	}
	return(threshnames)
}

mxExpectationNormal <-
  function(covariance, means = NA, dimnames = NA, thresholds = NA,
           threshnames = dimnames, ..., discrete = as.character(NA),
           discreteSpec = NULL) {
	prohibitDotdotdot(list(...))
	if (missing(covariance) || typeof(covariance) != "character") {
		stop("'covariance' argument is not a string (the name of the expected covariance matrix)")
	}
	if (!(single.na(means) || typeof(means) == "character")) {
		stop("Means argument is not a string (the name of the expected means matrix)")
	}
	if (is.na(means)) means <- as.integer(NA)
	if (single.na(thresholds)) thresholds <- as.character(NA)
	if (single.na(dimnames)) dimnames <- as.character(NA)
	if (!is.vector(dimnames) || typeof(dimnames) != 'character') {
		stop("'dimnames' argument is not a character vector")
	}
	if (length(thresholds) != 1) {
		stop("'thresholds' argument must be a single matrix or algebra name")
	}
	if (length(dimnames) == 0) {
		stop("'dimnames' argument cannot be an empty vector")
	}
	if (length(dimnames) > 1 && any(is.na(dimnames))) {
		stop("NA values are not allowed for 'dimnames' vector")
	}
	threshnames <- checkThreshnames(threshnames)
	return(new("MxExpectationNormal", covariance, means, dimnames, thresholds, threshnames,
             discrete, discreteSpec))
}

displayMxExpectationNormal <- function(expectation) {
	cat("MxExpectationNormal", omxQuotes(expectation@name), '\n')
	cat("$covariance :", omxQuotes(expectation@covariance), '\n')
	cat("$means :", omxQuotes(expectation@means), '\n')
	if (single.na(expectation@dims)) {
		cat("$dims : NA \n")
	} else {
		cat("$dims :", omxQuotes(expectation@dims), '\n')
	}
	if (single.na(expectation@thresholds)) {
		cat("$thresholds : NA \n")
	} else {
		cat("$thresholds :", omxQuotes(expectation@thresholds), '\n')
	}
	if (is.null(expectation@discreteSpec)) {
		cat("$discreteSpec : NULL \n")
	} else {
		cat("$discreteSpec :", expectation@discreteSpec, '\n')
	}
	if (single.na(expectation@discrete)) {
		cat("$discrete : NA \n")
	} else {
		cat("$discrete :", omxQuotes(expectation@discrete), '\n')
	}
	if (single.na(expectation@threshnames)) {
		cat("$threshnames : NA \n")
	} else {
		cat("$threshnames :", omxQuotes(expectation@threshnames), '\n')
	}
	invisible(expectation)
}


setMethod("print", "MxExpectationNormal", function(x,...) {
	displayMxExpectationNormal(x)
})

setMethod("show", "MxExpectationNormal", function(object) {
	displayMxExpectationNormal(object)
})
