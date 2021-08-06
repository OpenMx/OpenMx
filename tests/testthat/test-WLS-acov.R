library(testthat)
require(OpenMx)
context("WLS acov")
data(Bollen)

got <- mxGenerateData(Bollen[, 1:8], nrows=10)
omxCheckEquals(nrow(got), 10)

#--------------------------------------
# Set up  model matrices

manvar <- names(Bollen[, 1:8])

lval <- matrix(
	c(1, 0,
		1, 0,
		1, 0,
		1, 0,
		0, 1,
		0, 1,
		0, 1,
		0, 1),
	byrow=TRUE,
	ncol=2, nrow=8)
lfre <- matrix(as.logical(lval), ncol=2)
lfre[1, 1] <- FALSE
lfre[5, 2] <- FALSE
llab <- matrix(c(paste("lam", 1:4, sep=""), rep(NA, 8), paste("lam", 1:4, sep="")), ncol=2)


lx <- mxMatrix(name="Lam", values=lval, free=lfre, ncol=2, nrow=8, labels=llab, dimnames=list(manvar, c("F1", "F2")))


td <- mxMatrix(name="Theta", type="Symm", ncol=8,
							 values=
							 	c(.8,  0,  0,  0, .2,  0,  0,  0,
							 		.8,  0, .2,  0, .2,  0,  0,
							 		.8,  0,  0,  0, .2,  0,
							 		.8,  0,  0,  0, .2,
							 		.8,  0,  0,  0,
							 		.8,  0, .2,
							 		.8,  0,
							 		.8),
							 free=c(T,F,F,F,T,F,F,F,
							 			 T,F,T,F,T,F,F,
							 			 T,F,F,F,T,F,
							 			 T,F,F,F,T,
							 			 T,F,F,F,
							 			 T,F,T,
							 			 T,F,
							 			 T),
							 dimnames=list(manvar, manvar)
)
diag(td$labels) <- paste("var", 1:8, sep="")
selMat <- matrix(
	c(5,1,
		4,2,
		6,2,
		7,3,
		8,4,
		8,6), ncol=2, byrow=TRUE)
td$labels[selMat] <- paste("cov", c(51, 42, 62, 73, 84, 86), sep="")
td$labels[selMat[,2:1]] <- paste("cov", c(51, 42, 62, 73, 84, 86), sep="")

ph <- mxMatrix(name="Phi", type="Symm", ncol=2, free=T, values=c(.8, .2, .8), labels=paste("phi", c(1, 12, 2), sep=""), dimnames=list(c("F1", "F2"), c("F1", "F2")))


#--------------------------------------
# Set-up WLS model

wlsMod <- mxModel(
	"Test case for WLS Objective function from Bollen 1989",
	lx, ph, td,
	mxExpectationLISREL(LX=lx$name, PH=ph$name, TD=td$name),
	mxFitFunctionWLS(),
	mxData(Bollen[, 1:8], 'raw')
)

wlsRun <- mxRun(wlsMod)
omxCheckTrue(is.null(wlsRun$output$calculatedHessian))

dwlsMod <- mxModel(wlsMod, mxFitFunctionWLS("DWLS"))
dwlsRun <- mxRun(dwlsMod)

mxdw <- omxAugmentDataWithWLSSummary(mxd=wlsMod$data, type="DWLS")

dwlsMod2 <- dwlsMod
dwlsMod2$data <- mxData(
	mxdw$observedStats$cov, numObs = 75, means = NA, type = "acov", 
	#acov=diag(diag(mxdw$observedStats$acov)),
	fullWeight=mxdw$observedStats$asymCov * 75,
	acov=mxdw$observedStats$useWeight)
dwlsRun2 <- mxRun(dwlsMod2)
expect_equivalent(coef(dwlsRun) - coef(dwlsRun2),
                  rep(0, length(coef(wlsRun))))

dwlsMod3 <- dwlsMod
dwlsMod3$data <- mxData(
  observedStats = mxdw$observedStats, numObs = 75)
dwlsRun3 <- mxRun(dwlsMod3)
expect_equivalent(coef(dwlsRun) - coef(dwlsRun3),
                  rep(0, length(coef(wlsRun))))

mxw <- omxAugmentDataWithWLSSummary(mxd=wlsMod$data)

wlsMod2 <- wlsMod
wlsMod2$data <- mxData(
  mxw$observedStats$cov, numObs = 75, means = NA, type = "acov", 
  #acov=diag(diag(mxdw$observedStats$acov)),
  fullWeight=mxw$observedStats$asymCov * 75,
  acov=mxw$observedStats$useWeight)
wlsRun2 <- mxRun(wlsMod2)
expect_equivalent(coef(wlsRun) - coef(wlsRun2),
                  rep(0, length(coef(wlsRun))))

wlsMod3 <- wlsMod
wlsMod3$data <- mxData(
  observedStats = mxw$observedStats, numObs = 75)
wlsRun3 <- mxRun(wlsMod3)
expect_equivalent(coef(wlsRun) - coef(wlsRun3),
                  rep(0, length(coef(wlsRun))))
