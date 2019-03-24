 ===========================
 = check mxSE isnt broken =
 ===========================

# Created: 2019-03-24 04:57PM
# path = "~/bin/OpenMx/inst/models/failing/check_mxSE_works.R"
# Was Broken as of Version: 2.12.2.233 [GIT v2.12.2-233-ga7a310a]
# NOTE: When fixed, this check should be merged into passing/mxSE_test.R

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

# Try and get an SE for the parameter by expression, and then by label

# By expression "works", but should not warn that model has been changed...)
mxSE(X, m1) #  0.4024916

# Warning message:
# In assertModelRunAndFresh(object) :
#   MxModel 'one_is_the_loneliest_number' was modified since it was run.


# By label errors, and should not...

mxSE("X", m1)
# Treating first argument as character named entity in the model
# Error in convertForBackend(flatModel@compute, flatModel, model) :
#   Can only apply MxComputeJacobian to MxAlgebra *or* MxExpectation not "X"


