library(OpenMx)
library(mvtnorm)
library(testthat)
context("LatentAR-190605")

if (1) {
  # Not sure if this simulation exactly matches the model. Parameter
  # recovery might be better with a better simulation.

  set.seed(1)
  tsData <- NULL
  
  N <- 500
  muI <- runif(1)
  muS <- runif(1)
  muB <- 1       # 1.0 is no effect
  Ve <- .2
  VI <- .5
  VS <- .5
  IScov <- -.25
  VB <- .2
  for (n in 1:N) {
    got <- rmvnorm(1, c(muI,muS), sigma=matrix(c(VI,IScov,IScov,VS),2,2))
    I1 <- got[1]
    S1 <- got[2]
    B1 <- rnorm(1, muB, VB)
    part <- rep(0, 6)
    part[1] <- rnorm(1,0,.2) + I1
    for (x in 1:5) {
      prev <- part[x] * B1
      part[x+1] <- part[x+1] + prev + rnorm(1, 0, Ve)
    }
    part <- part + S1 * 0:5
    tsData <- rbind(tsData, part)
  }
  colnames(tsData) <- paste0('t',1:ncol(tsData))
  tsData <- as.data.frame(tsData)
}

# Create the vectors of variable names
tManifests <- colnames(tsData)
tLGCLatents <- c("I", "S")
tARLatent <- c("B")
tOps <- paste("op", c(1:5), sep="")

# This fits a latent growth curve to "prewhiten the time series" 
testARModel <- mxModel(model="testAR", type="RAM",
    manifestVars=tManifests,
    latentVars=c(tLGCLatents,tARLatent,tOps), #
    mxPath(from="I", to='t1', arrows=1, free=FALSE, values=1),
    mxPath(from="S", to=tManifests, arrows=1, free=FALSE, values=c(0:5)),
    mxPath(from="I", to="S", arrows=2, free=TRUE, values=0, labels="IScov"),
    mxPath(from="B", to=tOps, arrows=1, free=FALSE, values=1),
    mxPath(from=tManifests[1:5], to=tOps, arrows=1, free=FALSE, values=1),
    mxPath(from=tOps, to=tManifests[2:6], arrows=1, free=FALSE, values=1),
    
    mxPath(from=tManifests, arrows=2, free=TRUE, values=1, labels="Ve"),
    mxPath(from=c(tLGCLatents), arrows=2, free=TRUE, values=c(.53,.52), labels=c("VI", "VS")),
    mxPath(from=c(tARLatent), arrows=2, free=TRUE, values=.11, labels=c("VB")),

    mxPath(from="one", to=c(tLGCLatents), arrows=1, free=TRUE, values=c(.5,.4), labels=c("muI", "muS")),
    mxPath(from="one", to=c(tARLatent), arrows=1, free=TRUE, values=1, labels=c("muB")),
    mxPath(from="one", to=tOps, arrows=1, free=FALSE, values=1),
    mxData(observed=cov(tsData), type='cov', means=colMeans(tsData),
           numObs=nrow(tsData))
)

#print(colnames(testARModel$A))

testARModel$expectation$isProductNode <- colnames(testARModel$A) %in% tOps
#testARModel$expectation$isProductNode <- rep(FALSE, nrow(testARModel$A))

testARModelFit <- mxRun(testARModel)
summary(testARModelFit)

expect_equal(testARModelFit$output$fit, -1660.948, .01)
expect_equivalent(coef(testARModelFit)['muI'], muI, tolerance=.06)
expect_equivalent(coef(testARModelFit)['muS'], muS, tolerance=.1)
expect_equivalent(coef(testARModelFit)['muB'], muB, tolerance=.6)
expect_equivalent(coef(testARModelFit)['Ve'], Ve, tolerance=.2)
expect_equivalent(coef(testARModelFit)['VI'], VI, .4)
expect_equivalent(coef(testARModelFit)['VS'], VS, .2)
expect_equivalent(coef(testARModelFit)['IScov'], IScov, .04)
expect_equivalent(coef(testARModelFit)['VB'], VB, .15)

if (0) {
  library(ggplot2)
  library(reshape2)
  tsData$row <- 1:nrow(tsData)
  df <- melt(tsData, id.vars="row")
  ggplot(df) +
    geom_line(aes(x=variable, y=value, group=row), alpha=.1)
}
