library(OpenMx)

uMod <- mxModel(
	"upper", type="RAM",
	latentVars = c("eta"),
	mxData(type="raw", observed=data.frame(id=1L), primaryKey = 'id'),
	mxPath('eta', arrows=2, values=1))

for (numLower in 2:10) {
	lMod <- mxModel(
		"lower", type="RAM", uMod,
		manifestVars = 'y1',
		mxData(type="raw", observed=data.frame(y1=rnorm(numLower), id=1L)),
		mxPath('one', 'y1'),
		mxPath('y1', arrows=2, values=1),
		mxPath('upper.eta', 'y1', values=.1, joinKey='id'),
		mxComputeSequence(list(
			mxComputeOnce('fitfunction', 'fit'),
			mxComputeReportExpectation())))

	lMod <- mxRun(lMod, silent=TRUE)

	target <- cbind(matrix(1, numLower, 1), contr.helmert(numLower)[numLower:1,(numLower-1):1])
	rot <- -qr.Q(qr(target))

	ed <- lMod$expectation$debug
	
	# Only the sign of the first term matters.
	omxCheckCloseEnough(abs(c(lMod$data$observed$y1 %*% rot)),
			    abs(c(ed$g2$dataVec, ed$g1$dataVec)), 1e-6)
}
