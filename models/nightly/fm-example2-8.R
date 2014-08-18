# flexMIRT 2.0 example 2.8
library(OpenMx)
library(rpf)

ex2.8 <- expand.grid(Match=ordered(0:3), Identify=ordered(0:3))
ex2.8$freq <- c(71, 34, 30, 38, 30, 13, 15, 15, 13, 4, 15, 19, 43, 30, 25,  197)

spec <- list()
spec[1:2] <- rpf.nrm(outcomes=4, T.a="id", T.c="id")

ip.mat <- mxMatrix(name="item", values=c(1,1,rep(.1,5)), free=TRUE, nrow=7, ncol=2,
                   dimnames=list(names(rpf.rparam(spec[[1]])), colnames(ex2.8)[c(2,1)]))
ip.mat$free['alf1',] <- FALSE
ip.mat$values['alf2', 'Identify'] <- 0
ip.mat$free['alf2', 'Identify'] <- FALSE
ip.mat$labels[c('alf2','alf3'), 'Match'] <- 'eq1'
ip.mat$labels[c('gam1','gam2'), 'Match'] <- 'eq2'

imat <- ip.mat
imat$values['a',] <- 0
imat$free['a',] <- FALSE
imat$free[c('alf2', 'alf3'),] <- FALSE

m1 <- mxModel(model="ex28", ip.mat,
              mxData(observed=ex2.8, type="raw", sort = FALSE, numObs = sum(ex2.8$freq)),
              mxExpectationBA81(spec, weightColumn="freq"),
              mxFitFunctionML())

m1 <- mxRun(mxModel(m1, mxComputeSequence(list(
  mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson()),
  mxComputeOnce('fitfunction', 'information', 'meat'),
  mxComputeHessianQuality()))), silent=TRUE)

omxCheckCloseEnough(m1$output$fit, 2766.688, .01)

m1Sum <- summary(m1)
omxCheckCloseEnough(m1Sum$observedStatistics, 15, .1)
omxCheckCloseEnough(m1Sum$informationCriteria['AIC:','par'], 2784.69, .01)
omxCheckCloseEnough(summary(m1)$informationCriteria['BIC:','par'], 2824.14, .01)

fitstat <- m1$compute$steps[[1]]$output
omxCheckCloseEnough(m1$output$evaluations, 218, 5)
omxCheckCloseEnough(fitstat$EMcycles, 16, 1)
omxCheckCloseEnough(fitstat$totalMstep, 106, 5)

refModels <- mxNullModels(m1, run=TRUE)
nullm1 <- refModels$Independence
omxCheckCloseEnough(nullm1$output$fit, 2885.665, .01)
omxCheckCloseEnough(summary(nullm1)$informationCriteria['AIC:','par'], 2895.67, .01)
omxCheckCloseEnough(summary(nullm1)$informationCriteria['BIC:','par'], 2917.58, .01)
