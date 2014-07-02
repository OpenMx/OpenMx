# flexMIRT 2.0 example 2.8
library(OpenMx)
library(rpf)

ex2.8 <- expand.grid(Match=ordered(0:3), Identify=ordered(0:3))
ex2.8$freq <- c(71, 34, 30, 38, 30, 13, 15, 15, 13, 4, 15, 19, 43, 30, 25,  197)

spec <- list()
spec[1:2] <- rpf.nrm(outcomes=4, T.a="id", T.c="id")

ip.mat <- mxMatrix(name="item", values=c(1,1,rep(.1,5)), free=TRUE, nrow=7, ncol=2,
                   dimnames=list(names(rpf.rparam(spec[[1]])), colnames(ex2.8)[c(2,1)]))
ip.mat$values['alf2', 'Identify'] <- 0
ip.mat$free['alf2', 'Identify'] <- FALSE
ip.mat$labels[c('alf2','alf3'), 'Match'] <- 'eq1'
ip.mat$labels[c('gam1','gam2'), 'Match'] <- 'eq2'

imat <- ip.mat
imat$values['a',] <- 0
imat$free['a',] <- FALSE

nullm1 <- mxModel(model="ex28", imat,
              mxData(observed=ex2.8, type="raw", sort = FALSE),
              mxExpectationBA81(spec, weightColumn="freq"),
              mxFitFunctionML())
nullm1 <- mxRun(mxModel(nullm1, mxComputeEM('expectation', 'scores',
                                            mxComputeNewtonRaphson())), silent=TRUE)
omxCheckCloseEnough(nullm1$output$fit, 2885.665, .01)

m1 <- mxModel(model="ex28", ip.mat,
              mxData(observed=ex2.8, type="raw", sort = FALSE),
              mxExpectationBA81(spec, weightColumn="freq"),
              mxFitFunctionML())

m1 <- mxRun(mxModel(m1, mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson())), silent=TRUE)
omxCheckCloseEnough(m1$output$fit, 2766.688, .01)
fitstat <- m1$compute$output
omxCheckCloseEnough(m1$output$evaluations, 1616, 50)
omxCheckCloseEnough(fitstat$EMcycles, 16, 1)
omxCheckCloseEnough(fitstat$totalMstep, 297, 10)
