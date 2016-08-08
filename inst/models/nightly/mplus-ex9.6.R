# MPLUS: TWO-LEVEL CFA WITH CONTINUOUS FACTOR INDICATORS AND COVARIATES
# See https://www.statmodel.com/usersguide/chapter9.shtml

library(OpenMx)

set.seed(1)
ex96 <- suppressWarnings(try(read.table("models/nightly/data/ex9.6.dat")))
if (is(ex96, "try-error")) ex96 <- read.table("data/ex9.6.dat")

ex96$V8 <- as.integer(ex96$V8)
bData <- ex96[!duplicated(ex96$V8), c('V7', 'V8')]
colnames(bData) <- c('w', 'clusterID')
wData <- ex96[,-match(c('V7'), colnames(ex96))]
colnames(wData) <- c(paste0('y', 1:4), paste0('x', 1:2), 'clusterID')

bModel <- mxModel(
    'between', type="RAM",
    mxData(type="raw", observed=bData, primaryKey="clusterID"),
    latentVars = c("lw", "fb"),
    mxPath("one", "lw", labels="data.w", free=FALSE),
    mxPath("fb", arrows=2, labels="psiB"),
    mxPath("lw", 'fb', labels="phi1"))

wModel <- mxModel(
    'within', type="RAM", bModel,
    mxData(type="raw", observed=wData),
    manifestVars = paste0('y', 1:4),
    latentVars = c('fw', paste0("xe", 1:2)),
    mxPath("one", paste0('y', 1:4), values=runif(4),
	   labels=paste0("gam0", 1:4)),
    mxPath("one", paste0('xe', 1:2),
	   labels=paste0('data.x',1:2), free=FALSE),
    mxPath(paste0('xe', 1:2), "fw",
	   labels=paste0('gam', 1:2, '1')),
    mxPath('fw', arrows=2, values=1.1, labels="varFW"),
    mxPath('fw', paste0('y', 1:4), free=c(FALSE, rep(TRUE, 3)),
	   values=c(1,runif(3)), labels=paste0("loadW", 1:4)),
    mxPath('between.fb', paste0('y', 1:4), values=c(1,runif(3)),
	   free=c(FALSE, rep(TRUE, 3)), labels=paste0("loadB", 1:4),
	   joinKey="clusterID"),
    mxPath(paste0('y', 1:4), arrows=2, values=rlnorm(4),
	   labels=paste0("thetaW", 1:4)))

mle <- structure(c(
    0.9989, 0.9948, 1.0171, 0.9809, 0.9475, 1.0699,
    1.0139, 0.9799, -0.0829, -0.0771, -0.0449, -0.0299, 0.9728, 0.5105,
    0.9595, 0.9238, 0.9489, 0.361, 0.3445),
		 .Names = c("loadW2", "loadW3", "loadW4", "thetaW1",
		     "thetaW2", "thetaW3", "thetaW4", "varFW",
		     "gam01", "gam02", "gam03", "gam04", "gam11", "gam21",
		     "loadB2", "loadB3", "loadB4", "psiB", "phi1"))

if (1) {
	pt1 <- omxSetParameters(wModel, labels=names(mle), values=mle)
#	pt1$expectation$.forceSingleGroup <- TRUE
#	pt1$expectation$.rampart <- 0L
	plan <- mxComputeSequence(list(
	    mxComputeOnce('fitfunction', 'fit'),
#	    mxComputeNumericDeriv(checkGradient=FALSE,
#				  hessian=FALSE, iterations=2),
	    mxComputeReportDeriv(),
	    mxComputeReportExpectation()
	))
	pt1 <- mxRun(mxModel(pt1, plan))
	omxCheckCloseEnough(pt1$output$fit, 13088.373, 1e-2)
}

if (1) {
#  wModel <- mxRun(mxModel(wModel, mxComputeGradientDescent(verbose=2L)))
  wModel <- mxRun(wModel)
  summary(wModel)

  omxCheckCloseEnough(wModel$output$fit, 13088.373, 1e-2)
  omxCheckCloseEnough(mle[names(coef(wModel))], coef(wModel), 1e-3)
  omxCheckCloseEnough(wModel$expectation$debug$rampartUsage, 890)
} else {
	options(width=120)
	plan <- mxComputeSequence(list(
	    mxComputeOnce('fitfunction', 'fit'),
	    mxComputeNumericDeriv(checkGradient=FALSE,
				  hessian=FALSE, iterations=2),
	    mxComputeReportDeriv(),
	    mxComputeReportExpectation()
	))

	wModel$expectation$.rampart <- 2L
#	wModel$expectation$scaleOverride <- c(6, 1)
	rotated <- mxRun(mxModel(wModel, plan))
	
	wModel$expectation$.rampart <- 0L
	square <- mxRun(mxModel(wModel, plan))

	ex <- rotated$expectation
	eo <- ex$output
	ed <- ex$debug
	print(ed$rampartUsage)
	print(abs(rotated$output$fit - square$output$fit))
	print(max(abs(rotated$output$gradient - square$output$gradient)))
}
