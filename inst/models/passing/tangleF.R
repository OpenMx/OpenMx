library(OpenMx)

# NPSOL trips up on WLS for some reason
if (mxOption(NULL,"Default optimizer") == 'NPSOL') stop("SKIP")

mxOption(NULL, "Standard Errors", "No")

numManifestsPerLatent <- 5

mkModel <- function(shuffle, fellner) {
  varStruct <- expand.grid(l=1:3, i=1:numManifestsPerLatent)
  manifestVars <- paste0('i', apply(varStruct, 1, paste0, collapse=''))
  if (shuffle) {
    manifestVars <- sample(manifestVars, length(manifestVars))
  }
  latentVars <- paste0('l', unique(varStruct$l))
  allVars <- c(manifestVars, latentVars)
  if (shuffle) {
    allVars <- sample(allVars, length(allVars))
  }
  
  Fval <- diag(length(manifestVars))[,match(allVars, manifestVars)]
  Fval[is.na(Fval)] <- 0
  
  ta1 <- mxModel(
    model="tangle",
    mxMatrix("Full", length(manifestVars), length(allVars),
             values=Fval,
             dimnames=list(manifestVars, allVars), name="F"),
    mxMatrix("Symm", length(allVars), length(allVars),
             values=diag(length(allVars)),
             free=diag(length(allVars)) == 1,
             dimnames=list(allVars, allVars), name="S"),
    mxMatrix("Full", length(allVars), length(allVars),
             values=0,
             dimnames=list(allVars, allVars), name="A"),
    mxMatrix("Full", 1, length(allVars),
             free=!is.na(match(allVars, manifestVars)),
             dimnames=list(NULL, allVars), name="M"),
    mxExpectationRAM(M="M"),
    mxFitFunctionML(fellner=fellner))
  
  for (lx in 1:length(latentVars)) {
    lvar <- paste0('l', lx)
    ivar <- paste0(paste0('i', lx), 1:numManifestsPerLatent)
    ta1$A$values[ivar, lvar] <- 1
    ta1$A$free[ivar, lvar[-1]] <- TRUE
  }
  
  ta1$S$free[latentVars, latentVars] <- TRUE
  ta1$expectation$.maxDebugGroups <- 10L
  ta1
}

set.seed(1)
ta1 <- mxGenerateData(mkModel(FALSE, FALSE), nrow=25, returnModel=TRUE)

for (useFellner in c(TRUE,FALSE)) {
	fit1 <- mxRun(mxModel(mkModel(TRUE, useFellner), ta1$data))
	fit2 <- mxRun(mxModel(mkModel(FALSE, useFellner), ta1$data))

	omxCheckCloseEnough(fit1$output$fit - fit2$output$fit, 0, 1e-8)

	if (useFellner) {
		omxCheckEquals(names(fit1$expectation$debug$g1$fullMean[1:ncol(fit1$F)]),
			       paste0('tangle.', colnames(fit1$F)))

		omxCheckTrue(length(rle(fit1$expectation$debug$g1$latentFilter)$lengths) !=
				 length(rle(fit2$expectation$debug$g1$latentFilter)$lengths))
	}
}

if (0) {
  # F does two things:
  # drops out latent variables
  # permutes the manifest variables
  
  all(rownames(fit1$F$values) == fit1$expectation$.runDims)
  filter <- !apply(fit1$F$values, 2, function(x) all(x==0))
  
  map1 <- apply(fit1$F$values != 0, 2, function(x) ifelse(any(x), which(x), NA))
  man1 <- colnames(fit1$F)[filter]
  # permit F input order to F output order
  all(man1[ match(rownames(fit1$F), man1) ] == rownames(fit1$F))
  # permute data to F input order
  all(rownames(fit1$F)[ match(man1, rownames(fit1$F)) ] == man1)
  all(match(man1, rownames(fit1$F)) == map1[!is.na(map1)])
}

# ----------------------------------------

data("jointdata", package ="OpenMx", verbose= TRUE)
jointData <- jointdata

jointData[,c(2,4,5)] <- mxFactor(jointData[,c(2,4,5)], 
				 levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

mkModel <- function(shuffle, wls) {
  set.seed(shuffle)
	myData <- jointData
	if (shuffle) {
	  dperm <- sample.int(ncol(myData), ncol(myData))
		myData <- myData[, dperm]
	}

	thresh <- mxMatrix("Full", 3, 3, FALSE, 0, name="Th")

	thresh$free[,1] <- c(TRUE, FALSE, FALSE)
	thresh$values[,1] <- c(0, NA, NA)
	thresh$labels[,1] <- c("z2t1", NA, NA)

	thresh$free[,2] <- TRUE
	thresh$values[,2] <- c(-1, 0, 1)
	thresh$labels[,2] <- c("z4t1", "z4t2", "z4t3")

	thresh$free[,3] <- c(TRUE, TRUE, FALSE)
	thresh$values[,3] <- c(-1, 1, NA)
	thresh$labels[,3] <- c("z5t1", "z5t2", NA)

	colnames(thresh) <- paste0('z',c(2,4,5))

	manifestVars <- colnames(jointData)
	if (shuffle) manifestVars <- sample(manifestVars, length(manifestVars))

	latentVars <- 'l1'
	allVars <- c(manifestVars, latentVars)
	if (shuffle) allVars <- sample(allVars, length(allVars))
  
	Fval <- diag(length(manifestVars))[,match(allVars, manifestVars)]
	Fval[is.na(Fval)] <- 0
  
	freeMean <- !sapply(myData, is.factor)[match(allVars,colnames(myData))]
	freeMean <- !is.na(freeMean) & freeMean

	ta1 <- mxModel(
		model="tangle", thresh,
		mxData(myData, 'raw'),
		mxMatrix("Full", length(manifestVars), length(allVars),
			 values=Fval,
			 dimnames=list(manifestVars, allVars), name="F"),
		mxMatrix("Symm", length(allVars), length(allVars),
			 values=diag(length(allVars)),
			 free=diag(length(allVars)) == 1, lbound=0, ubound=5,
			 dimnames=list(allVars, allVars), name="S"),
		mxMatrix("Full", length(allVars), length(allVars),
			 values=0,
			 dimnames=list(allVars, allVars), name="A"),
		mxMatrix("Full", 1, length(allVars),
			 free=freeMean,
			 dimnames=list(NULL, allVars), name="M"),
		mxExpectationRAM(M="M", thresholds="Th"),
		mxComputeGradientDescent())
#		mxComputeOnce('fitfunction', 'fit'))
  
	lvar <- 'l1'
	ivar <- manifestVars
	ta1$A$values[ivar, lvar] <- 1
	ta1$A$lbound[ivar, lvar] <- 0
	ta1$A$free[ivar, lvar] <- TRUE
	ta1$S$free[lvar,lvar] <- FALSE
	ta1$M$values[1,'z1'] <- c(.1)

	if (wls) {
		ta1 <- mxModel(ta1, mxFitFunctionWLS('WLS'))
	} else {
		ta1 <- mxModel(ta1, mxFitFunctionML(jointConditionOn = "continuous"))
	}

	ta1
}

for (wls in c(FALSE)) {
	fit1 <- mxRun(mkModel(0, wls))  # MLE=2683.071 when wls=false
	fit2 <- mxRun(mkModel(1, wls))
	omxCheckCloseEnough(fit1$output$fit - fit2$output$fit, 0, 1e-6)
	fit3 <- mxRun(mkModel(2, wls))
	omxCheckCloseEnough(fit1$output$fit - fit3$output$fit, 0, 1e-6)
}

for (sx in 1:4) {
  for (wls in c(TRUE)) {
    fit1 <- mxRun(mkModel(0, wls))  # MLE=2683.071 when wls=false
    fit2 <- mxRun(mxModel(mkModel(sx, wls), fit1$data))
    fit3 <- mxRun(mkModel(sx, wls))
    fit4 <- mxRun(mxModel(mkModel(0, wls), fit3$data))
#    print(c(fit1$output$fit, fit2$output$fit, fit3$output$fit, fit4$output$fit))
    omxCheckCloseEnough(fit1$output$fit - fit2$output$fit, 0, 1e-6)
    omxCheckCloseEnough(fit1$output$fit - fit3$output$fit, 0, 1e-6)
    omxCheckCloseEnough(fit1$output$fit - fit4$output$fit, 0, 1e-6)
  }
}
