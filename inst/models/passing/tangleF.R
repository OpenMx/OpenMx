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

jointData <- suppressWarnings(try(read.table("models/passing/data/jointdata.txt", header=TRUE), silent=TRUE))
jointData <- read.table("data/jointdata.txt", header=TRUE)

jointData[,c(2,4,5)] <- mxFactor(jointData[,c(2,4,5)], 
				 levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

acovPerm <- function(wd, perm) {
  sz <- length(perm)
  tcount <- colSums(!is.na(wd$thresholds))
  names(tcount) <- NULL
  tstart <- cumsum(c(0,tcount))
  
  tperm <- rep(NA, sum(tcount))
  to <- 1
  thresholdColumns <- match(colnames(wd$thresholds), colnames(wd$means)) #oldCol
  newOrder <- order(perm[thresholdColumns]) # newOrder
  
  for (t1 in 1:length(newOrder)) {
    oldIndex <- newOrder[t1]
    for (cx in 1:tcount[oldIndex]) {
      tperm[to] <- tstart[oldIndex] + cx
      to <- to + 1
    }
  }
#  print(tperm)

  iperm <- rep(NA,length(perm))
  for (xx in 1:length(iperm)) iperm[perm[xx]] = xx
  
  part1 <- (sz * (sz+1))/2
  
  mm <- matrix(NA,5,5)
  mm[lower.tri(mm, diag = TRUE)] <- 1:part1
  mm[upper.tri(mm)] <- t(mm)[upper.tri(mm)]
  #  print(mm[perm,perm]-1)
  
  aperm <- 1:nrow(wd$acov)
  aperm[1:part1] <- vech(mm[iperm,iperm])
  aperm[(part1 + 1):(part1 + sz)] <-
    aperm[(part1 + 1):(part1 + sz)][iperm]
  aperm[(part1 + sz + 1):(part1 + sz + sum(tcount))] <-
    aperm[(part1 + sz + 1):(part1 + sz + sum(tcount))][tperm]
  aperm
}

if (1) {  # demonstrate the effect of data column permutation
  wlsData <- mxDataWLS(jointData, type='WLS')
  perm <- sample.int(5,5)
  wlsData2 <- mxDataWLS(jointData[,perm], type='WLS')

  aperm <- acovPerm(wlsData2, perm)

  omxCheckCloseEnough(max(abs(diag(wlsData2$acov)[aperm] - diag(wlsData$acov))), 0, 1e-6)
  omxCheckCloseEnough(max(abs(wlsData2$acov[aperm,aperm] - wlsData$acov)), 0, 1e-6)
}

mkModel <- function(shuffle, wls) {
	myData <- jointData
	if (shuffle) {
#	  dperm <- c( 4,   0,   2,   3,   1) + 1
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
	  md <- suppressWarnings(mxDataWLS(myData, type='WLS'))
		ta1 <- mxModel(ta1, md, mxFitFunctionWLS())
	} else {
		ta1 <- mxModel(ta1, mxData(myData, type="raw"), mxFitFunctionML(jointConditionOn = "continuous"))
	}

	ta1
}

for (wls in c(FALSE,TRUE)) {
	fit1 <- mxRun(mkModel(FALSE, wls))  # MLE=2683.071 when wls=false
	fit2 <- mxRun(mkModel(TRUE, wls))
	fit3 <- mxRun(mkModel(TRUE, wls))

	omxCheckCloseEnough(fit1$output$fit - fit2$output$fit, 0, 1e-6)
	omxCheckCloseEnough(fit1$output$fit - fit3$output$fit, 0, 1e-6)
}
