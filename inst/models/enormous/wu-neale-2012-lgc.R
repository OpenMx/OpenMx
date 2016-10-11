library(OpenMx)

mxOption(NULL, "Default optimizer", "SLSQP")

if (mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")

#print(mxOption(NULL, "Default optimizer"))

resVars      <- mxPath( from=c("x1","x2","x3","x4","x5"), arrows=2,
                        free=TRUE,  values = 1,
                        labels=c("residual","residual","residual","residual","residual") )
latVars      <- mxPath( from=c("intercept","slope"), arrows=2, connect="unique.pairs",
                       free=c(TRUE,FALSE,TRUE), values=c(1,0,1), labels=c("vari","cov","vars"))
intLoads     <- mxPath( from="intercept", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=1 )
sloLoads     <- mxPath( from="slope", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=seq(-2,2) )
manMeans     <- mxPath( from="one", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=0)
latMeans     <- mxPath( from="one", to=c("intercept", "slope"), arrows=1,
                        free=FALSE, values=0, labels=c("meani","means") )

growthCurveModel <- mxModel("Linear Growth Curve Model Path Specification",
                             type="RAM",
                             manifestVars=c("x1","x2","x3","x4","x5"),
                             latentVars=c("intercept","slope"),
                             resVars, latVars, intLoads, sloLoads,
                             manMeans, latMeans)

result <- expand.grid(rep=1:25000, adj=c(TRUE,FALSE), trueSvar=c(0,.3,.6),
                      interval=c(.95),
                      lbound=NA, val=NA, ubound=NA, retries=NA)
if (0) {
  load("/tmp/lgc-sim.rda")
}
bounds <- c('lbound','ubound')

for (rx in 1:nrow(result)) {
#for (rx in which(is.na(result$lbound))) {
  set.seed(result[rx,'rep'])

  true.svar <- result[rx,'trueSvar']
  ci.adj <- result[rx,'adj']

  growthCurveModel$S$values['slope','slope'] <- true.svar
  dset <- mxGenerateData(growthCurveModel, nrows = 150)
    m1 <- mxModel(growthCurveModel,
                  mxData(cov(dset), 'cov', colMeans(dset), nrow(dset)))

    m1$S$values['slope','slope'] <- .5
    if (ci.adj) {
      m1$S$lbound['slope','slope'] <- 0
    } else {
      m1$S$lbound['slope','slope'] <- NA
    }
    m1 <- mxModel(m1, mxCI('vars', boundAdj = ci.adj, interval=result[rx,'interval']))
    plan <- mxComputeSequence(list(
      GD=mxComputeGradientDescent(),
      CI=mxComputeConfidenceInterval(
        fitfunction="fitfunction",
        constraintType='ineq',
        plan=mxComputeTryHard(
          maxRetries=50L, scale=.05,
          plan=mxComputeGradientDescent(nudgeZeroStarts = FALSE,
                                        maxMajorIter = 150L)))))
    m1 <- mxModel(m1, plan)
    m1 <- try(mxRun(m1, intervals=TRUE, suppressWarnings=TRUE, silent=TRUE))
    if (is(m1, "try-error")) {
      print(paste("optimizer failed on", rx)) #52128
      next
    }

    detail <- m1$compute$steps[['CI']]$output$detail
    ci <- m1$output$confidenceIntervals
    result[rx,bounds] <- ci[1,bounds]
    result[rx,'val'] <- ci[1,'estimate']
    result[rx,'retries'] <- m1$compute$steps[['CI']]$plan$debug$retries
    
    if (rx %% 1000 == 0) {
      print(rx)
      save(result, file="/tmp/lgc-sim.rda")
    }
}

result$region <- NA
result[result[,'lbound'] <= result[,'trueSvar'] &
         result[,'trueSvar'] <= result[,'ubound'],'region'] <- 'M'
result[result[,'lbound'] > result[,'trueSvar'], 'region'] <- 'L'
result[result[,'ubound'] < result[,'trueSvar'], 'region'] <- 'U'

library(plyr)
resultSummary <- ddply(result, .(adj, trueSvar, interval), function(sim) {
  c(M=sum(sim$region == 'M'), L=sum(sim$region== 'L'), U=sum(sim$region=='U')) /
    nrow(sim)
})
print(resultSummary)

# how to estimate the Monte Carlo standard error?

if(0) {
  save(resultSummary, file="~/vcu/ci/lgc-sim.rda")
}
if(0) {
  lgc.sim <- result
  save(lgc.sim, file="~/vcu/ci/lgc-sim.rda")
}
