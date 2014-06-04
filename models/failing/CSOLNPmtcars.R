# models/failing/CSOLNPmtcars.R
# Not taking big enough steps?

library(OpenMx)
library(testthat)
# ===================
# = test with NPSOL =
# ===================

mxOption(NULL, "Default optimizer", "NPSOL")

# ====================================
# = make a simple independence model =
# ====================================
manifests = c("mpg", "disp", "gear")
m1 <- mxModel("ind", type = "RAM",
	manifestVars = manifests,
	mxPath(from = manifests, arrows = 2),
	mxPath(from = "one", to = manifests),
	mxData(mtcars[,manifests], type="raw")
)

# ======================================
# = get parameter estimated parameters =
# ======================================
m1 = mxRun(m1)
obtainedVars = summary(m1)$parameters[1:3,5]
expectedVars = as.numeric(diag(cov(mtcars[,manifests])))
mlexpected = expectedVars*(31/32) # n-1/n
omxCheckWithinPercentError(obtainedVars[1], mlexpected[1], percent = 0.1)
omxCheckWithinPercentError(obtainedVars[2], mlexpected[2], percent = 0.1)
omxCheckWithinPercentError(obtainedVars[3], mlexpected[3], percent = 0.1)

# ====================
# = switch to CSOLNP =
# ====================

mxOption(NULL, "Default optimizer", "CSOLNP")
manifests = c("mpg", "disp", "gear")
# re-build to set to crappy starts that are fine in NPSOL
m1 <- mxModel("ind", type = "RAM",
	manifestVars = manifests,
	mxPath(from = manifests, arrows = 2),
	mxPath(from = "one", to = manifests),
	mxData(mtcars[,manifests], type="raw")
)
m1 = mxRun(m1)
obtainedVars = summary(m1)$parameters[1:3,5]
expectedVars = as.numeric(diag(cov(mtcars[,manifests])))
mlexpected = expectedVars*(31/32) # n-1/n

testthat::expect_equal(obtainedVars[1], expected = mlexpected[1], tolerance = 0.01)
testthat::expect_equal(obtainedVars[2], expected = mlexpected[2], tolerance = 0.01)
testthat::expect_equal(obtainedVars[3], expected = mlexpected[3], tolerance = 0.01)
# e.g
# Error: obtainedVars[3] not equal to expectedVars[3]
# Mean relative difference: 34520.75