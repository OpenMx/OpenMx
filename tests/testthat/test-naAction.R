library(OpenMx)
library(testthat)
context("naAction")
suppressWarnings(RNGversion("3.5"))

data(demoOneFactor)
dof <- demoOneFactor

dof$x5 <- as.integer(dof$x5)  # test autoconversion to numeric

mask <- matrix(as.logical(rbinom(prod(dim(dof)), size = 1, .1)),
               nrow=nrow(dof), ncol=ncol(dof))
dof[mask] <- NA
manifests <- names(dof)
latents <- c("G")
model <- mxModel("OneFactor", 
                 type="LISREL",
                 manifestVars=list(exo=manifests), 
                 latentVars=list(exo=latents),
                 mxPath(from=latents, to=manifests),
                 mxPath(from=manifests, arrows=2),
                 mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
                 mxPath(from='one', to=manifests),
                 mxData(observed=dof, type="raw", naAction = 'fail'))
expect_error(mxRun(model),
             "contains at least one NA")

model$data$observed <- dof[rowSums(mask) == 0,]
fit <- mxRun(model)
fit1 <- fit$output$fit

model$data$observed <- dof
model$data$naAction <- 'omit'
d1 <- mxRun(model)
expect_equal(d1$output$fit, fit1, 1e-9)

set.seed(1)
b1 <- mxBootstrap(d1, 10)
expect_equal(sd(b1$compute$output$raw[,'fit']), 28, 1)

model$data$naAction <- 'exclude'
d2 <- mxRun(model)
expect_equal(d2$output$fit, fit1, 1e-9)

set.seed(1)
b2 <- mxBootstrap(d2, 10)
expect_equal(b1$compute$output$raw[,'fit'], b2$compute$output$raw[,'fit'], 1e-9)

# ----

numSets <- 4
ob <- list()
for (rep in 1:numSets) {
  v <- rnorm(nrow(dof))
  v[as.logical(rbinom(nrow(dof), size = 1, .1))] <- NA
  ob[[rep]] <- v
}
ob <- as.data.frame(ob, col.names=paste0('V',1:numSets))

m1 <- mxModel(
  model,
  mxFitFunctionWLS(allContinuousMethod = "marginals"),
  mxComputeLoop(list(
    LD=mxComputeLoadData("OneFactor", column="x5", method="data.frame",
                         byrow=FALSE, observed=ob),
    mxComputeSetOriginalStarts(),
    mxComputeOnce('fitfunction','fit'),
    CP=mxComputeCheckpoint(toReturn = TRUE, parameters=FALSE))))

m1 <- mxRun(m1)

m2 <- m1
m2$data$naAction <- 'omit'
m2 <- mxRun(m2)
expect_equal(m2$compute$steps$CP$log$objective, m1$compute$steps$CP$log$objective)

# ----

m3 <- m1
m3$data$observed[['freq']] <- 1L + rpois(nrow(dof), .5)
m3$data$frequency <- 'freq'
m3$data$naAction <- 'exclude'
m4 <- m3
m3 <- mxRun(m3)
m4$data$naAction <- 'omit'
m4 <- mxRun(m4)
expect_equal(m3$compute$steps$CP$log$objective,
             m4$compute$steps$CP$log$objective)
expect_equal(nrow(m3$data$observed), 500)
expect_equal(nrow(m4$data$observed), 300,10)

# ----

numSets <- 4
ob <- list()
for (rep in 1:numSets) {
  ob[[rep]] <- rnorm(nrow(dof))
}
ob <- as.data.frame(ob, col.names=paste0('V',1:numSets))

m1 <- mxModel(
  model,
  mxFitFunctionWLS(allContinuousMethod = "marginals"),
  mxComputeLoop(list(
    LD=mxComputeLoadData("OneFactor", column="x5", method="data.frame",
                         byrow=FALSE, observed=ob),
    mxComputeSetOriginalStarts(),
    mxComputeOnce('fitfunction','fit'),
    CP=mxComputeCheckpoint(toReturn = TRUE, parameters=FALSE))))

m1 <- mxRun(m1)

m2 <- m1
m2$data$naAction <- 'omit'
m2 <- mxRun(m2)
expect_equal(m2$compute$steps$CP$log$objective, m1$compute$steps$CP$log$objective)

# ----
