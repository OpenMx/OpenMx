library(OpenMx)

#mxOption(NULL, "Default optimizer", "SLSQP")
#print(mxOption(NULL, "Default optimizer"))

if (mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")

resVars      <- mxPath( from=c("x1","x2","x3","x4","x5"), arrows=2,
                        free=TRUE,  values = c(1,1,1,1,1),
                        labels=c("residual","residual","residual","residual","residual") )
latVars      <- mxPath( from=c("intercept","slope"), arrows=2, connect="unique.pairs",
                       free=c(TRUE,TRUE,FALSE), values=c(1,0,1),
                       labels=c("vari","cov","SlopeVarAlg[1,1]"))
intLoads     <- mxPath( from="intercept", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=c(1,1,1,1,1) )
sloLoads     <- mxPath( from="slope", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=seq(-2,2) )
manMeans     <- mxPath( from="one", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=runif(5,-.5,.5))
latMeans     <- mxPath( from="one", to=c("intercept", "slope"), arrows=1,
                        free=TRUE, values=runif(2, -.5,.5), labels=c("meani","means") )

SlopeVar <- mxMatrix("Full", 1,1, free=TRUE, values=0, name="SlopeVar")
SlopeVarAlg <- mxAlgebra(-SlopeVar, name="SlopeVarAlg")

growthCurveModel <- mxModel("sym",
                             type="RAM",
                             manifestVars=c("x1","x2","x3","x4","x5"),
                             latentVars=c("intercept","slope"),
                             resVars, latVars, intLoads, sloLoads,
                             manMeans, latMeans, SlopeVar, SlopeVarAlg)

result <- expand.grid(sign=c(1,-1), val=seq(-.1, .1, 0.005),
                     lbound=NA, ubound=NA, retries=NA)
bounds <- c('lbound','ubound')

cov1 <- mxGetExpected(growthCurveModel, "covariance")
mvec <- as.vector(mxGetExpected(growthCurveModel, "means"))
names(mvec) <- paste0('x',1:5)
growthCurveModel <- mxModel(growthCurveModel,
                            mxData(cov1, 'cov', mvec, 150))

interval <- c(.85,.95)

for (rx in 1:nrow(result)) {
  if (result[rx,'sign'] == 1) {
    SlopeVarAlg <- mxAlgebra(SlopeVar, name="SlopeVarAlg")
  } else {
    SlopeVarAlg <- mxAlgebra(-SlopeVar, name="SlopeVarAlg")
  }
  growthCurveModel <- mxModel(growthCurveModel, SlopeVarAlg)
  
    growthCurveModel$SlopeVar$values[1,1] <- result[rx,'sign'] * result[rx,'val']
    gcm <- mxRun(mxModel(growthCurveModel,
                         mxComputeSequence(list(
                           mxComputeOnce('expectation'),
                           mxComputeReportExpectation()))), silent=TRUE)

    cov1[,] <- gcm$expectation$output$covariance
    mvec[] <- gcm$expectation$output$mean
    m1 <- mxModel(growthCurveModel, mxData(cov1, 'cov', mvec, 150))
    
    if (any(eigen(cov1)$val < 0)) next

    if (result[rx,'sign'] == 1) {
      m1$SlopeVar$lbound[1,1] <- 0
    } else {
      m1$SlopeVar$ubound[1,1] <- 0
    }
    m1$SlopeVar$values[1,1] <- result[rx,'sign'] * 0.5
#    print(coef(m1))
    sinterval <- interval
    if (result[rx,'sign'] == -1) {
      sinterval <- rev(sinterval)
    }
    chkpt=FALSE
    m1 <- mxModel(m1, mxCI('SlopeVar', 'upper', interval=sinterval[1]))
    m2 <- mxModel(m1, mxCI('SlopeVar', 'lower', interval=sinterval[2]))
    m1 <- mxRename(m1, paste0('m',rx,'u'))
    m2 <- mxRename(m2, paste0('m',rx,'l'))
    m1 <- mxRun(m1, intervals=TRUE, suppressWarnings=TRUE, silent=TRUE,
                checkpoint=chkpt)
    m2 <- mxRun(m2, intervals=TRUE, suppressWarnings=TRUE, silent=TRUE,
                checkpoint=chkpt)
    
    detail <- m1$compute$steps[['CI']]$output$detail
    ci <- m1$output$confidenceIntervals
    result[rx,'ubound'] <- ci[1,'ubound']
    ci <- m2$output$confidenceIntervals
    result[rx,'lbound'] <- ci[1,'lbound']
    result[rx,'retries'] <- (m1$compute$steps[['CI']]$plan$debug$retries +
                               m2$compute$steps[['CI']]$plan$debug$retries)
}

if (FALSE) {
  library(ggplot2)
  library(gtable)
  library(grid)
  presult <- result
  presult$sign <- as.factor(presult$sign)
  p1 <- ggplot(presult) + geom_ribbon(aes(x=val, ymin=lbound, ymax=ubound, fill=sign), alpha=.3)
  p2 <- ggplot() + geom_point(data=presult[,c('val','retries','sign')], aes(x=val, y=retries, color=sign))
  pair <- rbind(ggplotGrob(p1), ggplotGrob(p2), size="first")
  grid.newpage()
  grid.draw(pair)
}

diff <- (result[result$sign == 1,c('lbound','ubound')] +
           result[result$sign == -1, c('ubound','lbound')])

omxCheckCloseEnough(sum(is.na(diff)), 0, 2)

omxCheckCloseEnough(max(abs(diff), na.rm=TRUE), 0, 1e-5)
