# ---------------------------------------------------------------------
# Program: OneFactorPathDemo.R
#  Author: Steve Boker
#    Date: Thu Jul 30 13:33:08 EDT 2009
#
# This program is the OpenMx one factor path model demo for the front page
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Thu Jul 30 13:33:11 EDT 2009
#      Created OneFactorPathDemo.R.
#
# ---------------------------------------------------------------------

require(OpenMx)
data(demoOneFactor)
manifestVars <- names(demoOneFactor)
factorModel <- mxModel("One Factor",
    mxMatrix("Full", 5, 1, values=0.2, free=T, name="A"),
    mxMatrix("Symm", 1, 1, values=1, free=F, name="L"),
    mxMatrix("Diag", 5, 5, values=1, free=T, name="U"),
    mxAlgebra(A %*% L %*% t(A) + U, name="R"),
    mxMLObjective("R", dimnames = manifestVars),
    mxData(cov(demoOneFactor), type="cov", numObs=500))
summary(mxRun(factorModel))

