# flexMIRT 2.0 example 2.8
library(OpenMx)
library(rpf)

ex2.8 <- expand.grid(Match=ordered(0:3), Identify=ordered(0:3))
ex2.8$freq <- c(71, 34, 30, 38, 30, 13, 15, 15, 13, 4, 15, 19, 43, 30, 25,  197)

spec <- list()
spec[1:2] <- list(rpf.nrm(outcomes=4, T.a="id", T.c="id"))

ip.mat <- mxMatrix(name="item", values=c(1,1,rep(.1,5)), free=TRUE, nrow=7, ncol=2,
                   dimnames=list(names(rpf.rparam(spec[[1]])), colnames(ex2.8)[c(2,1)]))
ip.mat$free['alf1',] <- FALSE
ip.mat$values['alf2', 'Identify'] <- 0
ip.mat$free['alf2', 'Identify'] <- FALSE
ip.mat$labels[c('alf2','alf3'), 'Match'] <- 'eq1'
ip.mat$labels[c('gam1','gam2'), 'Match'] <- 'eq2'

m1 <- mxModel(model="ex28", ip.mat,
              mxData(observed=ex2.8, type="raw", weight="freq"),
              mxExpectationBA81(spec),
              mxFitFunctionML())

m1 <- mxRun(mxModel(m1, mxComputeSequence(list(
  mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson()),
  mxComputeOnce('fitfunction', 'information', 'meat'),
  mxComputeHessianQuality()))), silent=TRUE)

omxCheckCloseEnough(m1$output$fit, 2766.688, .01)

m1Sum <- summary(m1)
omxCheckCloseEnough(m1Sum$observedStatistics, 15, .1)
omxCheckCloseEnough(m1Sum$informationCriteria['AIC:','par'], 2784.69, .01)
omxCheckCloseEnough(m1Sum$informationCriteria['BIC:','par'], 2791.64, .01)

fitstat <- m1$compute$steps[[1]]$output
omxCheckCloseEnough(m1$output$evaluations, 136, 5)
omxCheckCloseEnough(fitstat$EMcycles, 11, 1)
omxCheckCloseEnough(fitstat$totalMstep, 69, 5)

refModels <- mxRefModels(m1, run=TRUE)
omxCheckCloseEnough(refModels[[1]], c(2758.23, 0), .01)  #saturated
nullm1 <- refModels$Independence
omxCheckCloseEnough(nullm1$output$fit, 2885.665, .01)
omxCheckCloseEnough(summary(nullm1)$informationCriteria['AIC:','par'], 2895.67, .01)
omxCheckCloseEnough(summary(nullm1)$informationCriteria['BIC:','par'], 2899.52, .01)

# These fit indices are approximate because we use the number of
# statistics (30) for the df instead of the number of multinomial
# cells (32). This doesn't really matter though because the these
# statistics become inaccurate when the multinomial table is sparse
# and there is no known correction (Bock, Gibbons, & Muraki, 1998,
# p. 265).

m1Sum <- summary(m1, refModels=refModels)
omxCheckCloseEnough(m1Sum$CFI, .979, .01)
omxCheckCloseEnough(m1Sum$TLI, .965, .01)
omxCheckCloseEnough(m1Sum$RMSEA, .159, .01)
