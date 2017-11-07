# http://openmx.ssri.psu.edu/issue/2013/10/mxrename-fitted-model-objects-leaves-mxconstraint-objects-unrenamed

require(OpenMx)
 
aMat <- mxMatrix("Full", 1,1, name="aMat")
aCon <- mxConstraint(diag2vec(aMat)==0,name="aCon")
 
rTom  <- mxModel( "Tom", aMat, aCon)
summary(rTomRun <- mxRun(rTom))

rTomRun$aMat$values[1,1] <- .1
omxCheckWarning(logLik(rTomRun), "MxModel 'Tom' was modified since it was run.")

rNeal  <- mxRename(rTomRun, newname="Neal")
summary(rNeal)
