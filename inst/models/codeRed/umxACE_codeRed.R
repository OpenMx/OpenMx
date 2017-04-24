# All three optimizers give code RED
# SLSQP quits at start
# CSOLNP AND NPSOL appear to get good solutions

# ============================
# = How heritable is height? =
# ============================
# get a fresh build of umx
# requires umx 1.8 for inline optimizer switching
# install_github("tbates/umx")
require(umx)
oldJoint <- function(model){
	# https://github.com/OpenMx/OpenMx/commit/4fb4c0f190395f2b63b9710d986e547b835bfe92
	model$MZ$fitfunction$jointConditionOn <- 'old'
	model$DZ$fitfunction$jointConditionOn <- 'old'
	model = mxRun(model)
	umxSummary(model, std = FALSE)
	return(model)
}

data(twinData) # ?twinData from Australian twins.
# Pick the variables
selDVs = c("ht1", "ht2")
mzData <- twinData[twinData$zygosity %in% "MZFF", ]
dzData <- twinData[twinData$zygosity %in% "DZFF", ]

# -2ll used to be 9659
# All three optimizers code RED
m1 = umxACE(selDVs = selDVs, dzData = dzData, mzData = mzData, opt= "NPSOL")
# -2 × log(Likelihood) -11985.57 (df=4)
# |    |   a1|   c1|   e1|
# |:---|----:|----:|----:|
# |ht1 | 0.92| 0.14| 0.36|
# Warning message:
# In model 'ACE' Optimizer returned a non-zero status code 6. The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)
mold = oldJoint(m1) # Still code Red crazy estimates: |ht1 | 0.06| 0.01| 0.02|

m1 = umxACE(selDVs = selDVs, dzData = dzData, mzData = mzData, opt= "SLSQP")
# Code RED and start-value estimates
# 'log Lik.' -8963.089 (df=4)
# |    |   a1|   c1|   e1|
# |:---|----:|----:|----:|
# |ht1 | 0.58| 0.58| 0.58|
mold = oldJoint(m1) # Still code Red, crazy equal low estimates: |ht1 | 0.02| 0.02| 0.02|

m1 = umxACE(selDVs = selDVs, dzData = dzData, mzData = mzData, opt= "CSOLNP")
# Code red, but OK solution
# 'log Lik.' -11985.57 (df=4)
# |    |   a1|   c1|   e1|
# |:---|----:|----:|----:|
# |ht1 | 0.92| 0.14| 0.36|
mold = oldJoint(m1) # Still code Red, crazy equal low estimates: |ht1 | 0.06| 0.01| 0.02|
