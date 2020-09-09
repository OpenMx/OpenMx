library(testthat)
library(OpenMx)
context("Pearson selection")

t1 <- mxModel(
  "t1",
  mxMatrix("Symm", 3,3, FALSE, values=diag(3)+.2, name="c1"),
  mxMatrix("Symm", 3,3, free=c(F,T,rep(F,4)),
           values=diag(3)+.2, labels=paste0('p',1:6), name="c2"),
  mxAlgebra(mxPearsonSelCov(c1,c2), name="c3"),
  mxMatrix("Full", 3, 1, values=.1, name="m1"),
  mxAlgebra(mxPearsonSelMean(c1,c2,m1), name="m2"),
  mxExpectationNormal(covariance = "c3", means = "m2"),
  mxFitFunctionML())

trueSelCor <- -.2
t1 <- omxSetParameters(t1, "p2", values=trueSelCor)
if (0) {
  mxEval(c3, t1, compute = TRUE)
  mxEval(m2, t1, compute = TRUE)
}

set.seed(1)
dat <- mxGenerateData(t1, nrows = 800)

t2 <- mxModel(
  "t2", type="RAM",
  manifestVars = paste0('V',1:3),
  mxPath(paste0('V',1:3), arrows=2, values=1.2, free=FALSE),
  mxPath(paste0('V',1:3), arrows=2,
         connect = "unique.bivariate", values=.2, free=FALSE),
  mxPath('one',  paste0('V',1:3), values=.1, free=FALSE),
  mxPath('V1','V2', arrows=0, values=.2),
  mxData(dat, 'raw'))

t2 <- mxRun(t2)
summary(t2)

expect_equivalent(coef(t2), trueSelCor, .01)

expect_equivalent(t2$expectation$output$covariance,
                  mxGetExpected(t2, "covariance"), 1e-6)

expect_equivalent(t2$expectation$output$mean,
                  mxGetExpected(t2, "mean"), 1e-6)



# -- test multi-step TODO
