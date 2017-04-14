# ===========
# = HISTORY =
# ===========
# 2017-04-14 04:53PM Check by TBATES and still failing checks lines 48 and 49
# models/failing/CSOLNPmtcars.R
# CSOLNP Not taking big enough steps?

library(OpenMx)
# ===================
# = test with NPSOL =
# ===================

expectedVars = as.numeric(diag(cov(mtcars[,manifests])))
mlexpected = expectedVars*(31/32) # n-1/n

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
mxOption(NULL, "Default optimizer", "NPSOL")
m2 = mxRun(m1)
obtainedVars = summary(m2)$parameters[1:3,5]
omxCheckWithinPercentError(obtainedVars[1], mlexpected[1], percent = 0.1)
omxCheckWithinPercentError(obtainedVars[2], mlexpected[2], percent = 0.1)
omxCheckWithinPercentError(obtainedVars[3], mlexpected[3], percent = 0.1)

# ====================
# = switch to CSOLNP =
# ====================

mxOption(NULL, "Default optimizer", "CSOLNP")
m3 = mxRun(m1)
# In model 'ind' Optimizer returned a non-zero status code 5. The Hessian at the solution does not appear to be convex. See ?mxCheckIdentification for possible diagnosis (Mx status RED).
omxCheckTrue(mxCheckIdentification(m1)$status)

obtainedVars = summary(m3)$parameters[1:3,5]

omxCheckWithinPercentError(obtainedVars[1], mlexpected[1], percent = 0.01)
omxCheckWithinPercentError(obtainedVars[2], mlexpected[2], percent = 0.01)
omxCheckWithinPercentError(obtainedVars[3], mlexpected[3], percent = 0.01)
# e.g
# Error: obtainedVars[3] not equal to expectedVars[3]
# Mean relative difference: 34520.75