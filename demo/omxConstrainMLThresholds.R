# a simple one factor ordinal model
require(OpenMx)

data(myFADataRaw)

oneFactorOrd <- myFADataRaw[,c("z1", "z2", "z3")]
oneFactorOrd$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
oneFactorOrd$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
oneFactorOrd$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))

oneFactorModel <- mxModel("Common Factor Model Path Specification", 
	type="RAM",
	mxData(
		observed=oneFactorOrd,
		type="raw"
	),
	manifestVars=c("z1","z2","z3"),
	latentVars="F1",
	# residual variances
	mxPath(
		from=c("z1","z2","z3"),
		arrows=2,
		free=FALSE,
		values=c(1,1,1),
		labels=c("e1","e2","e3")
	),
	# latent variance
	mxPath(
		from="F1",
		arrows=2,
		free=TRUE,
		values=1,
		labels ="varF1"
	),
	# factor loadings
	mxPath(
		from="F1",
		to=c("z1","z2","z3"),
		arrows=1,
		free=c(FALSE,TRUE,TRUE),
		values=c(1,1,1),
		labels=c("l1","l2","l3")
	),
	# means
	mxPath(
		from="one",
		to=c("z1","z2","z3","F1"),
		arrows=1,
		free=FALSE,
		values=0,
		labels=c("meanz1","meanz2","meanz3","meanF")
	),
	# thresholds
	mxThreshold(vars=c("z1", "z2", "z3"),
		nThresh=c(1,1,2),
		free=TRUE,
		values=c(-1, 0, -.5, 1.2)
		)
) # close model

oneFactorCon <- omxConstrainMLThresholds(oneFactorModel)
oneFactorResults <- mxRun(oneFactorCon)
