library(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")

base <- mxModel(
	"OneFactorCov", type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests, values=rnorm(length(manifests))),
    mxPath(from=manifests, arrows=2, values=rlnorm(length(manifests))),
    mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
    mxPath(from = 'one', to = manifests, values=0),
	mxData(demoOneFactor, type="raw"))
base <- mxRun(base)

cmp1 <- mxRename(base, "cmp1")
cmp1$A$values['x1', 'G'] <- 0
cmp1$A$free['x1','G'] <- FALSE
cmp1 <- mxRun(cmp1)

cmp2 <- mxRename(base, "cmp2")
cmp2$M$values[,] <- 0
cmp2$M$free[1,'x1'] <- FALSE
cmp2 <- mxRun(cmp2)

cmp3 <- mxRename(cmp2, "cmp3")
cmp3$M$values[,] <- 0
cmp3$M$free[1,'x2'] <- FALSE
cmp3 <- mxRun(cmp3)

mxCompare(base, cmp1)  # p=0

omxCheckError(anova(cmp1),
	      "Compare model 'cmp1' with which other models?")

pgot <- anova(cmp1, cmp2, cmp3, base)
omxCheckEquals(is.na(pgot[,'p']), c(T,F,T,T))

set.seed(170623)
got <-  mxCompare(base, cmp1, boot=TRUE)
omxCheckCloseEnough(got[2,'p'], 0, 1/nrow(attr(got,'bootData')[[1]]))

pgot <- mxCompare(base, cmp2)

got <-  mxCompare(base, cmp2, replications=10)
omxCheckEquals(nrow(attr(got,'bootData')[[1]]), 10)

got <-  mxCompare(base, cmp2, previousRun = got, replications=500)
omxCheckCloseEnough(got[2,'p'], pgot[2,'p'], .01)

pgot <- mxCompareMatrix(list(base, cmp2, cmp3), 'p')
got <- mxCompareMatrix(list(base, cmp2, cmp3), 'p', replications=500)
omxCheckCloseEnough(cor(pgot[lower.tri(got)],
                        got[lower.tri(got)]), 1, .01)

