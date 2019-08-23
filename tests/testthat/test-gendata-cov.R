library(OpenMx)
library(testthat)
context("gendata-cov")

suppressWarnings(RNGversion("3.5"))
set.seed(1)

m1 = mxModel("testCovPower", type="RAM",
             mxData(observed = matrix(c(1, .3, .3, 1), nrow = 2, dimnames=list(c("X", "Y"), c("X", "Y"))), numObs = 200, type="cov"),
             manifestVars = c("X", "Y"),
             mxPath("X", to = "Y", value= .3),
             mxPath(c("X", "Y"), arrows=2, value= 1)
)

m2 = mxGenerateData(m1, nrows=200, returnModel=T)
omxCheckEquals(m2$data$type, "cov")
omxCheckEquals(m2$data$numObs, m1$data$numObs)
c1 <- m2$data$observed - mxGetExpected(m2, 'covariance')
omxCheckCloseEnough(range(abs(c1)), rep(.03, 2), .2)
omxCheckTrue(is.na(m2$data$means))

m3 = mxGenerateData(m1, nrows=200, returnModel=T, empirical=TRUE)
omxCheckEquals(m3$data$type, "cov")
omxCheckEquals(m3$data$numObs, m1$data$numObs)
omxCheckCloseEnough(m3$data$observed, mxGetExpected(m3, 'covariance'))
omxCheckTrue(is.na(m3$data$means))
