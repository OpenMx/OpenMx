library(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")

base <- mxModel(
	"OneFactorCov", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests, values=0,  free=FALSE, labels=paste0('l',1:length(manifests))),
    mxPath(from=manifests, arrows=2, values=rlnorm(length(manifests)), lbound=.01),
    mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
    mxPath(from = 'one', to = manifests, values=0, free=TRUE, labels=paste0('m',1:length(manifests))),
	mxData(demoOneFactor, type="raw"))
base <- mxRun(base)

set.seed(1)
got <- mxParametricBootstrap(base, paste0('l', 1:length(manifests)),
                             alternative="two.sided",
                             alpha=0.05,
                             correction="none",
                             replications=10)

omxCheckEquals(got['l2', 'note'], "< 1/10")

got2 <- mxParametricBootstrap(base, paste0('l', 1:length(manifests)), "two.sided",
                             replications=100,
                             previousRun=got)
omxCheckCloseEnough(got2[,'p'], c(.32,.32,.32,.25,.32), .12)

omxCheckEquals(attr(got,'bootData')[5,],
               attr(got2,'bootData')[5,])

got3 <- mxParametricBootstrap(base, paste0('l', 1:length(manifests)), "two.sided",
                              previousRun=got2)
omxCheckCloseEnough(got3[,'p'], c(.36,.36,.36,.325,.36), .08)

gc()
before <- proc.time()[['elapsed']]
got4 <- mxParametricBootstrap(base, paste0('l', 1:length(manifests)), "two.sided",
                              previousRun=got3)
elapsed <- proc.time()[['elapsed']] - before
omxCheckCloseEnough(elapsed, 0, .5)
