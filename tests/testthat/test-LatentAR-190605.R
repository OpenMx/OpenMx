library(OpenMx)
library(testthat)
context("LatentAR-190605")

if (0) {
  # This might be wrong TODO
  set.seed(1)
  tsData <- NULL
  
  N <- 500
  muI <- runif(1)
  print(muI)
  muS <- runif(1)
  print(muS)
  muB <- 1       # 1.0 is no effect
  for (n in 1:N) {
    I1 <- rnorm(1, muI, .2)
    S1 <- rnorm(1, muS, .2)
    B1 <- rnorm(1, muB, .2)
    part <- rep(0, 6)
    part[1] <- rnorm(1,0,.2)
    for (x in 1:5) {
      prev <- part[x] * B1
      #    prev <- 0
      part[x+1] <- part[x+1] + prev + rnorm(1, 0, .2)
    }
    part <- part + I1 * rep(1,6) + S1 * 0:5
    tsData <- rbind(tsData, part)
  }
  colnames(tsData) <- paste0('t',1:ncol(tsData))
}

data("latentAR_sim")
tsData <- latentAR_sim

# Create the vectors of variable names
tManifests <- colnames(tsData)
tLGCLatents <- c("I", "S")
tARLatent <- c("B")
tOps <- paste("op", c(1:5), sep="")

# This fits a latent growth curve to "prewhiten the time series" 
testARModel <- mxModel(model="testAR", type="RAM",
    manifestVars=tManifests,
    latentVars=c(tLGCLatents,tARLatent,tOps), #
    mxPath(from="I", to=tManifests, arrows=1, free=FALSE, values=1),
    mxPath(from="S", to=tManifests, arrows=1, free=FALSE, values=c(0:5)),
    mxPath(from="B", to=tOps, arrows=1, free=FALSE, values=1),
    mxPath(from=tManifests[1:5], to=tOps, arrows=1, free=FALSE, values=1),
    mxPath(from=tOps, to=tManifests[2:6], arrows=1, free=FALSE, values=1),
    
    mxPath(from=tManifests, arrows=2, free=TRUE, values=1, labels="Ve"),
    mxPath(from=c(tLGCLatents), arrows=2, free=TRUE, values=c(.93,.92), labels=c("VI", "VS")),
    mxPath(from=c(tARLatent), arrows=2, free=TRUE, values=.91, labels=c("VB")),

    mxPath(from="one", to=c(tLGCLatents), arrows=1, free=TRUE, values=c(.5,.4), labels=c("muI", "muS")),
    mxPath(from="one", to=c(tARLatent), arrows=1, free=TRUE, values=.9, labels=c("muB")),
    mxPath(from="one", to=tOps, arrows=1, free=FALSE, values=1),
    mxData(observed=cov(tsData), type='cov', means=colMeans(tsData),
           numObs=nrow(tsData))
)

#print(colnames(testARModel$A))

testARModel$expectation$isProductNode <- colnames(testARModel$A) %in% tOps
#testARModel$expectation$isProductNode <- rep(FALSE, nrow(testARModel$A))

testARModelFit <- mxRun(testARModel)
summary(testARModelFit)

#cat(deparse(round(coef(testARModelFit), 3)))
onyx <- c(Ve = 0.99, VI = 0.958, VS = 1.01, VB = 0.187, muI = 0.442,  muS = 0.487, muB = 0.915)
expect_equal(coef(testARModelFit), onyx, tolerance=1e-3)

if (0) {
  library(ggplot2)
  library(reshape2)
  tsData$row <- 1:nrow(tsData)
  df <- melt(tsData, id.vars="row")
  ggplot(df) +
    geom_line(aes(x=variable, y=value, group=row), alpha=.1)
}
