# MPLUS: TWO-LEVEL CFA WITH CONTINUOUS FACTOR INDICATORS AND COVARIATES

library(OpenMx)

set.seed(1)
ex96 <- suppressWarnings(try(read.table("models/nightly/data/ex9.6.dat")))
if (is(ex96, "try-error")) ex96 <- read.table("data/ex9.6.dat")

ex96$V8 <- as.integer(ex96$V8)
bData <- ex96[!duplicated(ex96$V8), c('V7', 'V8')]
colnames(bData) <- c('w', 'clusterID')
wData <- ex96[,-match(c('V7'), colnames(ex96))]
colnames(wData) <- c(paste0('y', 1:4), paste0('x', 1:2), 'clusterID')

bModel <- mxModel('between', type="RAM",
                  mxData(type="raw", observed=bData, primaryKey="clusterID"),
                  latentVars = c("lw", "fb"),
                  mxPath("one", "lw", labels="data.w", free=FALSE),
                  mxPath("fb", arrows=2, labels="v_fb"),
                  mxPath("lw", 'fb', labels="lw"))

wModel <- mxModel('within', type="RAM", bModel,
                  mxData(type="raw", observed=wData, sort=FALSE),  #[abs(wData$clusterID - 41)<= 25,]
                  manifestVars = paste0('y', 1:4),
                  latentVars = c('fw', paste0("xe", 1:2)),
                  mxPath("one", paste0('y', 1:4), values=runif(4), labels=paste0("ym", 1:4)),
                  mxPath("one", paste0('xe', 1:2), labels=paste0('data.x',1:2), free=FALSE),
                  mxPath(paste0('xe', 1:2), "fw", labels=paste0('xl', 1:2)),
                  mxPath('fw', arrows=2, values=1.1, labels="v_fw"),
                  mxPath('fw', paste0('y', 1:4), free=c(FALSE, rep(TRUE, 3)),
                         values=runif(3), labels=paste0("ly", 1:4)),
                  mxPath('between.fb', paste0('y', 1:4), values=runif(3),
                         free=c(FALSE, rep(TRUE, 3)), labels=paste0("lb", 1:4),
                         joinKey="clusterID"),
                  mxPath(paste0('y', 1:4), arrows=2, values=rlnorm(4), labels=paste0("v_y", 1:4)))

mle <- structure(c(0.999, 0.995, 1.017, 0.973, 0.51, 0.981, 0.948, 1.07,
		   1.014, 0.98, 0.96, 0.924, 0.949, -0.083, -0.077, -0.045, -0.03,  0.344, 0.361),
		 .Names = c("ly2", "ly3", "ly4", "xl1", "xl2",  "v_y1", "v_y2", "v_y3", "v_y4",
		     "v_fw", "lb2", "lb3", "lb4",  "ym1", "ym2", "ym3", "ym4", "lw", "v_fb"))

if (1) {
	pt1 <- omxSetParameters(wModel, labels=names(mle), values=mle)
	pt1$expectation$.forceSingleGroup <- TRUE
	pt1$expectation$rampart <- 0L
	pt1 <- mxRun(mxModel(pt1, mxComputeOnce('fitfunction', 'fit')))
	
	omxCheckCloseEnough(pt1$output$fit, 13088.373, 1e-2)
}

if (0) {
#  wModel <- mxRun(mxModel(wModel, mxComputeGradientDescent(verbose=2L)))
  wModel <- mxRun(wModel)
  summary(wModel)

  omxCheckCloseEnough(wModel$output$fit, 13088.373, 1e-2)
  omxCheckCloseEnough(mle, coef(wModel), 1e-2)
  omxCheckCloseEnough(wModel$expectation$debug$rampartUsage, 890)
} else {
	options(width=120)
	plan <- mxComputeSequence(list(
	    mxComputeOnce('fitfunction', 'fit'),
	    mxComputeNumericDeriv(checkGradient=FALSE, hessian=FALSE, iterations=2),
	    mxComputeReportDeriv(),
	    mxComputeReportExpectation()
	))

	wModel$expectation$rampart <- 2L
#	wModel$expectation$scaleOverride <- c(6, 1)
	rotated <- mxRun(mxModel(wModel, plan))
	
	wModel$expectation$rampart <- 0L
	square <- mxRun(mxModel(wModel, plan))

	ex <- rotated$expectation
	eo <- ex$output
	ed <- ex$debug
	print(ed$rampartUsage)
	print(abs(rotated$output$fit - square$output$fit))
	print(max(abs(rotated$output$gradient - square$output$gradient)))
}
