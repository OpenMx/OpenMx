library(OpenMx)

nContPerFactor <- 4
nOrdPerFactor <- 1
nVarPerFactor <- nContPerFactor + nOrdPerFactor
nFact <- 1
latents <- rawToChar(as.raw(as.integer(charToRaw("A")) + 1:nFact - 1), multiple=T)
manifests <- apply(expand.grid(prefix=latents, 1:(nContPerFactor + nOrdPerFactor)),
                   1, paste, collapse="")
ordinals <- apply(expand.grid(prefix=latents, 1:nOrdPerFactor),
                  1, paste, collapse="")
nGroups <- 3L

mkGroup <- function(name) {
  big <- mxModel(
    name, type="RAM",
    latentVars = latents,
    manifestVars = manifests,
    mxMatrix(nrow=2, ncol=nOrdPerFactor*nFact, values=1:2/3,
             free=FALSE, dimnames=list(NULL, ordinals), name="thresholds"))
  
  for (fx in latents) {
    big <- mxModel(
      big,
      mxPath("one", paste0(fx, 1:nVarPerFactor), values=rnorm(nVarPerFactor, sd = .2)),
      mxPath(paste0(fx, 1:nVarPerFactor), arrows=2, values=1, labels="err"),
      mxPath(fx, paste0(fx, 1:nVarPerFactor),
             values=runif(nVarPerFactor, .25,.5),
             labels=paste0("l",1:nVarPerFactor)),
      mxPath(fx, arrows=2, free=FALSE, values=1))
    if (fx == 'A') next
    big <- mxModel(
      big,
      mxPath('A', fx, values=runif(1,.25,.5)))
  }
  
  big$expectation$thresholds <- 'thresholds'
  big
}

container <- mxModel("mg", mxFitFunctionMultigroup(paste0('g', 1:nGroups)))
for (gx in 1:nGroups) {
  container <- mxModel(container,
                       mkGroup(paste0("g",gx)))
}
container <- omxAssignFirstParameters(container)
trueCoef <- coef(container)

set.seed(123)
container <- mxGenerateData(container, nrows=300, returnModel = TRUE)

ml <- mxModel(name="ml", container)
ml <- mxRun(ml)

omxCheckCloseEnough(max(abs(coef(ml) - trueCoef)), 0, .2)

wls <- mxModel(name="wls", container)

for (gx in 1:nGroups) {
  grp <- wls[[paste0("g",gx)]]
  grp <- mxModel(grp, mxDataWLS(grp$data$observed))
  wls <- mxModel(wls, grp)
}
wls <- mxRun(wls)

omxCheckCloseEnough(max(abs(coef(wls) - trueCoef)), 0, .2)

omxCheckCloseEnough(cor(coef(ml), coef(wls)), 1, .01)
