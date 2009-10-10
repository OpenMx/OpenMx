# -----------------------------------------------------------------------
# Program: OneFactorMatrixDemo.R  
#  Author: Steve Boker
#    Date: 08 01 2009 
#
# OpenMx one factor matrix model demo from front page of website
# 
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

data(demoOneFactor)

factorModel <- mxModel(name ="One Factor")
    mxMatrix(type="Full", nrow=5, ncol=1, free=T, values=0.2, name="A")
    mxMatrix(type="Symm", nrow=1, ncol=1, free=T, values=1, name="L")
    mxMatrix(type="Diag", nrow=5, ncol=5, free=T, values=1, name="U")
    mxAlgebra(expression=A %*% L %*% t(A) + U, name="R")
    mxMLObjective(covariance="R", dimnames=names(demoOneFactor))
    mxData(observed=cov(demoOneFactor), type="cov", numObs=500)

factorModelFit <- mxRun(factorModel)
summary(factorModelFit)
