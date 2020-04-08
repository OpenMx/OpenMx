library(OpenMx)
library(testthat)
context("omp")

skip_if_not(imxHasOpenMP())

Sys.setenv(OMP_NUM_THREADS = '1')
mxOption(key='Number of Threads', value=2)

mData = matrix (1)
dimnames(mData) = list(c("X"), c("X"))

m1 = mxModel("one_is_the_loneliest_number", type="RAM",
	manifestVars = "X",
	mxPath(from="X", to = "X", arrows=2, lbound=0, labels= "X"),
	mxData(mData, type="cov", numObs = 10)
)
expect_error(mxRun(m1), "2 threads requested.")
