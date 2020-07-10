#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2020-07-09 17:20:21
# Filename: MxPredict
# Purpose: Create an S3 predict method for MxModel objects
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

##' \code{predict} method for \code{MxModel} objects
##' 
##' @param object an MxModel object from which predictions are desired
##' @param newdata an optional \code{data.frame} object. See details.
##' @param interval character indicating what kind of intervals are desired.  'none' gives no intervals, 'confidence', gives confidence intervals, 'prediction' gives prediction intervals.
##' @param method character the method used to create the predictions.  See details.
##' @param level the confidence or predictions level, ignored if not using intervals
##' @param type character the type of thing you want predicted: latent variables or manifest variables.
##' @param ... further named arguments
##' 
##' @details
##' The \code{newdata} argument is either a \code{data.frame} or \code{MxData} object.  In the latter case is replaces the data in the top level model.  In the former case, it is passed as the \code{observed} argument of \code{mxData} with \code{type='raw'} and must accept the same further arguments as the data in the model passed in the \code{object} argument.
##' 
##' The available methods for prediction are 'ML', 'WeightedML', 'Regression', and 'Kalman'.  See the help page for \code{\link{mxFactorScores}} for details on the first three of these.  The 'Kalman' method uses the Kalman filter to create predictions for state space models.
predict.MxModel <- function(object,
                              newdata=NULL,
                              interval=c('none', 'confidence', 'prediction'),
                              method=c('ML', 'WeightedML', 'Regression', 'Kalman'),
                              level=.95,
                              type=c('latent', 'observed'),
                              ...
                             ){
	method <- match.arg(arg=tolower(method), choices=tolower(c('ML', 'WeightedML', 'Regression', 'Kalman')))
	interval <- match.arg(arg=tolower(interval), choices=c('none', 'confidence', 'prediction'))
	type <- match.arg(arg=tolower(type), choices=c('latent', 'observed'))
	# TODO handle nonnull newdata
	# TODO handle multigroup models intelligently and recursively
	if("MxExpectationStateSpace" %in% class(object$expectation)){
		if(method != 'kalman'){
			stop("State space model prediction is currently only implemented for method='Kalman'")
		}
		ks <- mxKalmanScores(object)
		est <- ks$xPredicted
		# TODO ensure that est is always a matrix
		se <- t(sqrt(matrix(apply(ks$PPredicted, 3, diag), nrow=dim(ks$PPredicted)[1], ncol=dim(ks$PPredicted)[3])))
	} else {
		stop("predict method is currently only implemented for simple state space models")
	}
	pred <- array(NA, dim=c(nrow(est), ncol(est), 3))
	pred[,,1] <- est
	pred[,,2] <- est - qnorm(1 - (1-level)/2)*se
	pred[,,3] <- est + qnorm(1 - (1-level)/2)*se
	return(pred)
}

#a <- predict.MxModel(srun, method='Kalman')
#plot(a[,1,1], type='l', ylim=c(-3, 3))
#lines(1:201, a[,1,2], lty=2)
#lines(1:201, a[,1,3], lty=2)


#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



#------------------------------------------------------------------------------

