# -----------------------------------------------------------------------
# Program: OneFactorMatrixDemo.R  
#  Author: Steve Boker
#    Date: 07 30 2009 
#
# OpenMx one factor matrix model demo for front page of website
# 
# Revision History
#   Hermine Maes -- 02 22 2010 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

data(demoOneFactor)
manifestVars <- names(demoOneFactor)

factorModel <- mxModel("One Factor",
    mxMatrix(type="Full", nrow=5, ncol=1, values=0.2, free=TRUE, name="A"),
    mxMatrix(type="Symm", nrow=1, ncol=1, values=1, free=FALSE, name="L"),
    mxMatrix(type="Diag", nrow=5, ncol=5, values=1, free=TRUE, name="U"),
    mxAlgebra(expression=A %*% L %*% t(A) + U, name="R"),
    mxMLObjective(covariance="R", dimnames=manifestVars),
    mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
)

factorFit <- mxRun(factorModel)
summary(factorFit)

