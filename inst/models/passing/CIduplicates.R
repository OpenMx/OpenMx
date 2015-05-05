library(OpenMx)
data(myFADataRaw, package="OpenMx")
manifests = paste0("x",1:3)
myFADataRaw = myFADataRaw[, manifests]
latents = c("G")
m1 <- mxModel("m1", type="RAM",
	manifestVars = manifests, latentVars   = latents,
	mxPath(from = latents, to = manifests),
	mxPath(from = manifests, arrows = 2, labels = paste0(manifests, "_resid")),
	mxPath(from = latents, arrows = 2, free = F, values = 1), # latents fixed at 1
	mxData(cov(myFADataRaw, use="complete"), type = "cov", numObs = nrow(myFADataRaw))
)
m1$S$lbound <- .1
m1 = mxRun(m1)

tmp = mxRun(mxModel(m1, mxCI(c("x1_resid","S[1,1]"))), intervals = T)
omxCheckCloseEnough(nrow(tmp$output$confidenceIntervals), 1)

omxCheckCloseEnough(tmp$output$confidenceIntervals[1, c('lbound', 'ubound')],
                    c(.280, .430), .01)
