# MPLUS: TWO-LEVEL CFA WITH CONTINUOUS FACTOR INDICATORS AND COVARIATES

library(OpenMx)

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
                  mxData(type="raw", observed=wData, sort=FALSE),
                  manifestVars = paste0('y', 1:4),
                  latentVars = c('fw', paste0("xe", 1:2)),
                  mxPath("one", paste0('y', 1:4), values=1, labels=paste0("ym", 1:4)),
                  mxPath("one", paste0('xe', 1:2), labels=paste0('data.x',1:2), free=FALSE),
                  mxPath(paste0('xe', 1:2), "fw", labels=paste0('xl', 1:2)),
                  mxPath('fw', arrows=2, values=1, labels="v_fw"),
                  mxPath('fw', paste0('y', 1:4), free=c(FALSE, rep(TRUE, 3)),
                         values=1, labels=paste0("ly", 1:4)),
                  mxPath('between.fb', paste0('y', 1:4), values=1,
                         free=c(FALSE, rep(TRUE, 3)), labels=paste0("lb", 1:4),
                         joinOn="clusterID"),
                  mxPath(paste0('y', 1:4), arrows=2, values=1, labels=paste0("v_y", 1:4)))

if (0) {
  wModel <- omxSetParameters(
    wModel, labels=names(coef(wModel)),
    values=c(.999, .995, 1.017, .973, .51, .981, .947, 1.07, 1.014, .98,
             .96, .924, .949, -.083, -0.077, -.045, -.03, .344, .361))
}

if (1) {
#  wModel <- mxRun(mxModel(wModel, mxComputeGradientDescent(verbose=2L)))
  wModel <- mxRun(wModel)
  summary(wModel)
} else {
  #coef(wModel)
  wModel <- mxRun(mxModel(wModel,
                          mxComputeSequence(list(
                            mxComputeOnce('fitfunction','fit'),
                            mxComputeReportExpectation()))))
  wModel$output$fit / -2
  me <- wModel$expectation$output
  me$means[1:4]                 # model implied means of first row
  me$covariance[1:4,1:4]  # covariance of first row
}
