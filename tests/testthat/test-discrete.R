library(testthat)
context("discrete")
library(OpenMx)

# to test:
# Normal, RAM, LISREL
# ML / WLS
# with and without regular thresholds
# ordered factor vs raw count

factorModel <- mxModel(
  "One Factor",
  mxMatrix(nrow=1, ncol=5, free=FALSE, values=0, name="M"),
  mxMatrix("Full", 5, 1, values=0.8,
           free=TRUE, name="A"),
  mxMatrix("Symm", 1, 1, values=1,
           free=FALSE, name="L"),
  mxMatrix("Diag", 5, 5, values=1,
           free=TRUE, name="U"),
  mxAlgebra(A %*% L %*% t(A) + U, name="R"),
  mxMatrix(nrow=5, ncol=3,
           values=c(4, -2, 1, 1.1, NA,
                    5, -2, 1, .6, NA,
                    6, -2, 2, 4, .5),
           dimnames=list(c(), paste0('x',1:3)),
           name="D"),
  mxExpectationNormal(covariance = "R",
                      dimnames = paste0('x',1:5),
                      discrete = "D",
                      means = "M"),
  mxFitFunctionML())

mxGetExpected(factorModel, "thresholds")

factorModel <- mxGenerateData(factorModel, 400, returnModel = TRUE)

# convert to raw counts
# factorModel$data$observed$x1 <-
#   unclass(factorModel$data$observed$x1) - 1L

#factorModel$D$free[2,] <- TRUE
factorModel$D$free[4:5,] <- !is.na(factorModel$D$values[4:5,])

fit <- mxRun(factorModel)
summary(fit)
fit$D

q()

# ---------------

library(OpenMx)

manifests <- paste0('x',1:5)
latents <- c("G")
factorModel <- mxModel(
  "One Factor",
  type="RAM",
  manifestVars = manifests,
  latentVars = latents,
  mxPath(from=latents, to=manifests,values=0.8),
  mxPath(from=manifests, arrows=2,values=1),
  mxPath(from=latents, arrows=2,
         free=FALSE, values=1.0),
  mxMarginalPoisson(paste0("x",1:2), c(4,5), c(.6,.7)),
  mxMarginalNegativeBinomial("x3", 6, 4, .5))

round(mxGetExpected(factorModel, "thresholds"),3)

mxGenerateData(factorModel, 20)
