# ================================================================
# = check mxSE and other functions work with a 1-parameter model =
# ================================================================

# Created: 2019-03-24 04:57PM
# path = "~/bin/OpenMx/inst/models/passing/check_mxSE_works.R"

library(testthat)
context("check1_var_RAM_works")
library(OpenMx)

# ============================================
# = Make a super-simple 1 variance RAM model =
# ============================================
mData = matrix (1)
dimnames(mData) = list(c("X"), c("X"))

m1 = mxModel("one_is_the_loneliest_number", type="RAM",
	manifestVars = "X",
	mxPath(from="X", to = "X", arrows=2, lbound=0, labels= "X"),
	mxData(mData, type="cov", numObs = 10)
)
m1 = mxRun(m1)

# SE by expression with no warnings or errors
omxCheckCloseEnough(mxSE(X, m1), 0.4024916, .01)
omxCheckWarning(mxSE(X, m1), message=NA)

# SE by label with no warnings or errors
omxCheckCloseEnough(mxSE(X, m1), mxSE("X", m1), 1e-4)
foo <- "X"
omxCheckCloseEnough(mxSE(X, m1), mxSE(foo, m1, forceName=T), 1e-4)
omxCheckWarning(mxSE("X", m1), message=NA)

# Confint with 1 or more parameters in the model
expect_message(omxCheckCloseEnough(confint(m1), c(.111,1.688), 1e-3),
               "Wald type confidence")
