library(testthat)
library(OpenMx)
context("Pearson selection")

skip_if(.Platform$OS.type=="windows" && .Platform$r_arch=="i386")

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

expect_error(mxModel(t2, mxPath('one', 'V1', arrows=0)),
             "means path must be a single-headed arrow")

t2 <- mxRun(t2)
summary(t2)

expect_equivalent(coef(t2), trueSelCor, .01)

expect_equivalent(t2$expectation$output$covariance,
                  mxGetExpected(t2, "covariance"), 1e-6)

expect_equivalent(t2$expectation$output$mean,
                  mxGetExpected(t2, "mean"), 1e-6)

# ----------------------------------

t3 <- mxModel(
  "t3",
  mxMatrix("Symm", 4,4, FALSE, values=diag(4)+.2, name="c1"),
  mxMatrix("Symm", 4,4, free=c(F,T,rep(F,8)),
           values=diag(4)+.2, labels=c(NA,'p1',rep(NA,8)), name="c2"),
  mxAlgebra(mxPearsonSelCov(c1,c2), name="c3"),
  mxMatrix("Full", 4, 1, values=.1, name="m1"),
  mxAlgebra(mxPearsonSelMean(c1,c2,m1), name="m2"),
  mxExpectationNormal(covariance = "c3", means = "m2"),
  mxFitFunctionML())

trueSelCor <- c(-.2, -.3)
t3 <- omxSetParameters(t3, "p1", values=trueSelCor[1])

t3 <- mxModel(t3,
              mxMatrix("Symm", 4,4,free=c(F,F,T,rep(F,7)),
                       values=mxEval(c3, t3, compute = TRUE),
                       labels=c(NA,NA,'p2',rep(NA,7)), name="c4"),
              mxAlgebra(mxPearsonSelCov(c3,c4), name="c5"),
              mxAlgebra(mxPearsonSelMean(c3,c4,m2), name="m3"),
              mxExpectationNormal(covariance = "c5", means = "m3"))

t3 <- omxSetParameters(t3, "p2", values=trueSelCor[2])
#mxEval(m3, t3, compute = TRUE)

set.seed(1)
dat <- mxGenerateData(t3, nrows = 800)

t4 <- mxModel(
  "t4", type="RAM",
  manifestVars = paste0('V',1:4),
  mxPath(paste0('V',1:4), arrows=2, values=1.2, free=FALSE),
  mxPath(paste0('V',1:4), arrows=2,
         connect = "unique.bivariate", values=.2, free=FALSE),
  mxPath('one',  paste0('V',1:4), values=.1, free=FALSE),
  mxPath('V1','V2', arrows=0, values=.2),
  mxPath('V1','V3', arrows=0, values=.2, step = 2),
  mxData(dat, 'raw'))

t4 <- mxRun(t4)
summary(t4)

expect_equivalent(coef(t4), trueSelCor, .03)

expect_equivalent(t4$expectation$output$covariance,
                  mxGetExpected(t4, "covariance"), 1e-6)

expect_equivalent(t4$expectation$output$mean,
                  mxGetExpected(t4, "mean"), 1e-6)
