library(OpenMx)

if (mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")

#mxOption(NULL, "Default optimizer", "SLSQP")
#print(mxOption(NULL, "Default optimizer"))

resVars      <- mxPath( from=c("x1","x2","x3","x4","x5"), arrows=2,
                        free=TRUE,  values = c(1,1,1,1,1),
                        labels=c("residual","residual","residual","residual","residual") )
latVars      <- mxPath( from=c("intercept","slope"), arrows=2, connect="unique.pairs",
                       free=TRUE, values=c(1,0,1), labels=c("vari","cov","vars"))
intLoads     <- mxPath( from="intercept", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=c(1,1,1,1,1) )
sloLoads     <- mxPath( from="slope", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=seq(-2,2) )
manMeans     <- mxPath( from="one", to=c("x1","x2","x3","x4","x5"), arrows=1,
                        free=FALSE, values=runif(5,-.5,.5))
latMeans     <- mxPath( from="one", to=c("intercept", "slope"), arrows=1,
                        free=TRUE, values=runif(2, -.5,.5), labels=c("meani","means") )

growthCurveModel <- mxModel("Linear Growth Curve Model Path Specification",
                             type="RAM",
                             manifestVars=c("x1","x2","x3","x4","x5"),
                             latentVars=c("intercept","slope"),
                             resVars, latVars, intLoads, sloLoads,
                             manMeans, latMeans)

result <- expand.grid(adj=c(TRUE,FALSE), val=seq(-.1, .1, 0.002),
                     lbound=NA, ubound=NA, retries=NA)
bounds <- c('lbound','ubound')

for (rx in 1:nrow(result)) {
#    rx=125
    growthCurveModel$S$values['slope','slope'] <- result[rx,'val']

    cov1 <- mxGetExpected(growthCurveModel, "covariance")
    if (any(eigen(cov1)$val < 0)) next
    mvec <- as.vector(mxGetExpected(growthCurveModel, "means"))
    names(mvec) <- paste0('x',1:5)
    m1 <- mxModel(growthCurveModel,
                  mxData(cov1, 'cov', mvec, 150))

    m1$S$values['slope','slope'] <- .5
    if (result[rx,'adj']) {
      m1$S$lbound['slope','slope'] <- 0
    } else {
      m1$S$lbound['slope','slope'] <- NA
    }
    m1 <- mxModel(m1, mxCI('vars', boundAdj = result[rx,'adj']))
    m1 <- mxRun(m1, intervals=TRUE, suppressWarnings=TRUE, silent=TRUE)

    detail <- m1$compute$steps[['CI']]$output$detail
    ci <- m1$output$confidenceIntervals
    result[rx,bounds] <- ci[1,bounds]
    result[rx,'retries'] <- m1$compute$steps[['CI']]$plan$debug$retries
}

omxCheckCloseEnough(table(is.na(result[!result$adj,'lbound']))[[1]],
               100,1)
omxCheckCloseEnough(sum(diff(result[!result$adj,'lbound']) > 0, na.rm = TRUE),
                    99, 1)
omxCheckCloseEnough(table(is.na(result[!result$adj,'ubound']))[[1]],
                    100,1)
omxCheckCloseEnough(sum(diff(result[!result$adj,'ubound']) > 0, na.rm = TRUE),
                    99, 1)
omxCheckEquals(fivenum(result[!result$adj, 'retries'])[c(1,3)],
               c(2,2))

omxCheckCloseEnough(table(is.na(result[result$adj,'lbound']))[[1]],
                    100,1)
omxCheckCloseEnough(sum(diff(result[result$adj,'lbound']) > 0, na.rm = TRUE),
                    38, 3)
omxCheckCloseEnough(table(is.na(result[result$adj,'ubound']))[[1]],
                    101,1)
omxCheckCloseEnough(sum(diff(result[result$adj,'ubound']) > 0, na.rm = TRUE),
                    98, 1)
omxCheckEquals(fivenum(result[result$adj, 'retries'])[c(1,3)],
               c(2,2))

if (0) {
  library(ggplot2)
  library(gtable)
  library(grid)
  p1 <- ggplot(result) + geom_ribbon(aes(x=val, ymin=lbound, ymax=ubound, fill=adj), alpha=.3) + ylim(-.02,.1)
  p2 <- ggplot() + geom_point(data=result[,c('val','retries','adj')], aes(x=val, y=retries, color=adj))
  pair <- rbind(ggplotGrob(p1), ggplotGrob(p2), size="first")
  grid.newpage()
  grid.draw(pair)
}
