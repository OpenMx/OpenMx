# ===============================
# = Model fails on being re-run =
# ===============================

library(OpenMx)
# =================
# = 1. Make data set =
# =================
set.seed(159)
xdat <- data.frame(a=rnorm(10, mean=4.2), b=1:10)

# =========================================
# = 2. Make a model with a row objective  =
# =========================================
cmod <- mxModel(name='Estimation Row Model with Missingness',
	mxData(observed=xdat, type='raw'),
	mxMatrix(values=.75, ncol=2, nrow=1, free=TRUE, name='M'),
	mxAlgebra(omxSelectCols(M, existenceVector), name='fM'),
	mxAlgebra((filteredDataRow-fM)%^%2, name='rowAlgebra'),
	mxAlgebra(sum(rowResults), name='reduceAlgebra'),
	mxFitFunctionRow(
		rowAlgebra='rowAlgebra',
		reduceAlgebra='reduceAlgebra',
		dimnames=c('a', 'b'))
)

# ==========
# = 3. Run =
# ==========
cmod <- mxRun(cmod)

# ==============
# = Now re-run =
# ==============
cmod <- mxRun(cmod)

mxVersion()
# OpenMx version: 2.7.9.58 [GIT v2.7.9-58-ga116bb8]
# R version: R version 3.3.3 (2017-03-06)
# Platform: x86_64-apple-darwin13.4.0
# MacOS: 10.12.5
# Default optimiser: CSOLNP
# Warning message:
# In model 'Estimation Row Model with Missingness' Optimizer returned a non-zero status code 6. The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)

# Back in 28 March 2015, this yielded an error
# Error: The filteredDataRow cannot have name 'Estimation Row Model with Missingness.filteredDataRow' because this entity already exists in the model

# Expected: runs quickly, returning something near the already-reached solution.
# Obtained: fails with error copied above