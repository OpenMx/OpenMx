require(OpenMx)
aMat <- mxMatrix("Full", 1,1, name="aMat")
aCon <- mxConstraint(diag2vec(aMat)==0,name="aCon")
 
mx101a13  <- mxModel( "mx101a13", aMat, aCon)
summary(mx101a13run <- mxRun(mx101a13))
mx101a15  <- mxRename(mx101a13run, newname="mx101a15")
summary(mx101a15)
