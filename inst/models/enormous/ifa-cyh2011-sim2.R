# This is a replication of Cai, Yang, & Hansen (2011) simulation study #2.

#options(error = utils::recover)
library(OpenMx)
library(rpf)
library(mvtnorm)

spec <- list()
spec[1:10] <- rpf.drm(factors=3)
spec[c(5,11:15)] <- rpf.nrm(factors=3, outcomes = 4, T.a=diag(3), T.c=diag(3))

correct <- matrix(NA, 9, length(spec))
correct[1,] <- c(1.37, 1.42, 1.57, 1.24, 1.32, 0.83, 1.74, 1.72, 1.23, 1.44,
                 2.04, 2.02, 1.35, 0.92, 1.84)
correct[2,1:5] <- c(1.4, 0.88, 0.59, 0.9, .88)
correct[3,11:15] <- .95
correct[4,c(1:4, 6:10)] <- c(-.27, -.79,-.92, -1.03, -1.27, -2.05, -2.36, -2.31, -2.79)
correct[5,c(1:4, 6:10)] <- -1.1   # approx logit(.25)
correct[4:5,5] <- c(2.49, 2.1)
correct[7:9,5] <- c(-.91, -3.4, -3.27)
correct[7:9,11] <- c(-3.87, -6.81, -8.04)
correct[7:9,12] <- c(-.11, -2.01, -2.58)
correct[7:9,13] <- c(-.9, -2.83,-4.68)
correct[7:9,14] <- c(-1.2,-4.34,-6.63)
correct[7:9,15] <- c(.11, -.05, -.06)
free <- !is.na(correct)

starting <- correct
starting[1,] <- 1
starting[2,1:5] <- 1
starting[3,11:15] <- 1
starting[4,c(1:4, 6:10)] <- 0
starting[5,c(1:4, 6:10)] <- round(logit(.1),2)  # should not start at MLE
starting[c(4:5,7:9), 5] <- 0
starting[7:9, 11:15] <- 0

setFixedParam <- function(mat) {
  mat[6,c(1:4, 6:10)] <- logit(1)
  mat[3,1:10] <- 0
  mat[2,6:15] <- 0
  mat[6, 5] <- 3
  mat[4:6, 11:15] <- 1:3
  mat
}

correct <- setFixedParam(correct)
starting <- setFixedParam(starting)

mkModel <- function(startValues, eq) {
  imat <- mxMatrix(name="item", values=startValues, free=free,
                   dimnames=list(names(rpf.rparam(spec[[15]])),
                                 paste("i",1:15,sep="")))
  glabels <- paste("g", 1:9, sep="")
  imat$labels[5,c(1:4, 6:10)] <- glabels
  
  if (eq) {
    imat$labels[3, 11:15] <- "eq1"
  }
  
  gaussPrior <- mxMatrix(name="gaussPrior", nrow=1, ncol=length(glabels),
                         free=TRUE, labels=glabels, values=startValues[5,1])
  
  gaussM <- mxMatrix(name="gaussM", nrow=1, ncol=length(glabels), values=-1.1)
  gaussSD <- mxMatrix(name="gaussSD", nrow=1, ncol=length(glabels), values=.5)
  
  gaussFit <- mxAlgebra(sum(log(2*pi) + 2*log(gaussSD) +
                              (gaussPrior-gaussM)^2/gaussSD^2), name="gaussFit")
  gaussGrad <- mxAlgebra(2*(gaussPrior - gaussM)/gaussSD^2, name="gaussGrad",
                         dimnames=list(c(),gaussPrior$labels))
  gaussHess <- mxAlgebra(vec2diag(2/gaussSD^2), name="gaussHess",
                         dimnames=list(gaussPrior$labels, gaussPrior$labels))
  
  gaussModel <- mxModel(model="gaussModel", gaussPrior, gaussM, gaussSD,
                        gaussFit, gaussGrad, gaussHess,
                        mxFitFunctionAlgebra("gaussFit", gradient="gaussGrad", hessian="gaussHess"))
  
  itemModel <- mxModel(model="itemModel", imat,
                       mxData(type="raw", observed=rpf.sample(3000, spec, correct)),
                       mxExpectationBA81(spec, qpoints=21, qwidth=5),
                       mxFitFunctionML())
  
  simModel <- mxModel(model="sim2", itemModel, gaussModel,
                      mxFitFunctionMultigroup(paste(c("itemModel", "gaussModel"),
                                                    "fitfunction", sep=".")),
                      mxComputeEM("itemModel.expectation", "scores",
                                  mxComputeNewtonRaphson()))
  simModel
}

seModel <- function() {
  m <- mkModel(starting, FALSE)
  m$compute$mstep$tolerance <- 1e-10  # for S-EM
  m
}

if (file.exists("models/enormous/lib/stderrlib.R")) {
  source("models/enormous/lib/stderrlib.R")
} else if (file.exists("lib/stderrlib.R")) {
  source("lib/stderrlib.R")
} else {
  stop("Cannot find stderrlib.R")
}

mcModel <- function() mkModel(correct, FALSE)  # speed up convergence

name <- "cyh2011-sim2"
getMCdata(name, mcModel, correct[free], maxCondNum=NA)

omxCheckCloseEnough(max(abs(mcBias)), .1258, .001)
omxCheckCloseEnough(norm(mcBias, "2"), .29957, .001)
omxCheckCloseEnough(log(det(mcHessian)), 318.3, .1)

detail <- testPhase(seModel, 500, methods=c('estepH', 'tian', 'agile', 'meat', 'oakes'))
asem <- studyASEM(seModel)
smooth <- checkSmoothness(seModel)

rda <- paste(name, "-result.rda", sep="")
save(detail, asem, smooth, file=rda)

