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

# Try and get an SE for the parameter by expression, and then by label

omxCheckWarning(mxSE(X, m1), NA) #  0.4024916

omxCheckCloseEnough(mxSE(X, m1), mxSE("X", m1), 1e-4)

foo <- "X"
omxCheckCloseEnough(mxSE(X, m1), mxSE(foo, m1, forceName=T), 1e-4)
