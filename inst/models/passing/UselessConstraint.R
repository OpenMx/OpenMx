# http://openmx.psyc.virginia.edu/issue/2014/05/memory-leak-when-running-ram-model-constraint

library(OpenMx)

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModelPath <- mxModel("OneFactorPath",
                           type="RAM",
                           manifestVars = manifests,
                           latentVars = latents,
                           mxPath(from=latents, to=manifests,
                                  labels=paste("l",1:5,sep="")),
                           mxPath(from=manifests, arrows=2),
                           mxPath(from=latents, arrows=2,
                                  free=FALSE, values=1.0),
                           mxData(cov(demoOneFactor), type="cov",
                                  numObs=500),
                           mxAlgebra(S[6,6],name="GV"),
                           mxConstraint(GV-1==0,name="pointless"),
                           mxConstraint(GV>0,name="morePointless"))
#factorModelPath <- mxOption(factorModelPath,"Checkpoint Directory","C:/Work/OpenMx_dev/")
#factorModelPath <- mxOption(factorModelPath,"Checkpoint Units","evaluations")
#factorModelPath <- mxOption(factorModelPath,"Checkpoint Count",1)
factorFit <- try(mxRun(factorModelPath), silent = TRUE)

if (is(factorFit, "try-error")) {
  # good
} else {
  # Any constraints that show up here by mistake will have a zero gradient.
  omxCheckTrue(all(factorFit$output$gradient != 0))
  omxCheckCloseEnough(sqrt(sum(factorFit$output$gradient^2)), 0, .01)
}
