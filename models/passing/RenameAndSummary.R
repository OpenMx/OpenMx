# http://openmx.psyc.virginia.edu/issue/2013/10/mxrename-fitted-model-objects-leaves-mxconstraint-objects-unrenamed

# This test exists because MxSummary uses the information in runstate.
# Discard this test as soon as MxSummary is rewritten to use the
# information in the model instead of runstate.

require(OpenMx)
 
aMat <- mxMatrix("Full", 1,1, name="aMat")
aCon <- mxConstraint(diag2vec(aMat)==0,name="aCon")
 
rTom  <- mxModel( "Tom", aMat, aCon)
summary(rTomRun <- mxRun(rTom))
rNeal  <- mxRename(rTomRun, newname="Neal")
summary(rNeal)
