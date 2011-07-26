## Test model for multiple latent layers
# Two latents, three manifests.  Both latents predict all three manifests, and one 
# latent also predicts the other latent
# Written by Daniel Hackett, June 2011

require(OpenMx)
dataTest <- rbind(c(4,2,2),c(2,4,2),c(2,2,4))
manifests <- c('X','Y','Z')
colnames(dataTest) <- manifests
rownames(dataTest) <- manifests
latents <- c('G', 'H')
factorModel <- mxModel("Two Factor Multilayer",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from='G', to=manifests,value=c(2,1,1),free=FALSE),
	  mxPath(from='H', to=manifests,value=c(1,1,1),free=FALSE),
	  mxPath(from='H', to='G', value=0.1, free=FALSE),
      mxPath(from=manifests, arrows=2,value=c(1,1,1), labels=c('E1','E1','E1')),
      mxPath(from=latents, arrows=2,
            values=1.0),
      mxData(dataTest, type="cov",
            numObs=100))
mxTestPPML(factorModel)