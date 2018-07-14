library(OpenMx)

if (mxOption(NULL,"Default optimizer") == 'NPSOL') stop("SKIP")

# get data
fp <- system.file("models/passing/data/jointdata.txt", package = "OpenMx", mustWork = TRUE)
jointData <- read.table(fp, header=TRUE)

# specify ordinal columns as ordered factors
jointData[,c(2,4,5)] <- mxFactor(jointData[,c(2,4,5)], 
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))
	
satCov <- mxMatrix("Symm", 5, 5,
	free=TRUE, values=diag(5), name="C")
satCov$free[2,2] <- FALSE
satCov$free[4,4] <- FALSE
satCov$free[5,5] <- FALSE

loadings <- mxMatrix("Full", 1, 5,
	free=TRUE, values=1, name="L", lbound=0)
loadings$ubound[1,4] <- 2
loadings$ubound[1,5] <- 2
	
resid <- mxMatrix("Diag", 5, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=.5, name="U")
	
means <- mxMatrix("Full", 1, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=0, name="M")
	
thresh <- mxMatrix("Full", 3, 3, FALSE, 0, name="T")

thresh$free[,1] <- c(TRUE, FALSE, FALSE)
thresh$values[,1] <- c(0, NA, NA)
thresh$labels[,1] <- c("z2t1", NA, NA)

thresh$free[,2] <- TRUE
thresh$values[,2] <- c(-1, 0, 1)
thresh$labels[,2] <- c("z4t1", "z4t2", "z4t3")

thresh$free[,3] <- c(TRUE, TRUE, FALSE)
thresh$values[,3] <- c(-1, 1, NA)
thresh$labels[,3] <- c("z5t1", "z5t2", NA)
	
omxCheckError(mxExpectationNormal("C", "M", dimnames=names(jointData),
                                  thresholds="T", threshnames=c("z2", "z4", "z2")),
              (paste("'threshnames' argument contains 'z2' more than once. \nIf you are having problems with Doppelgangers",
                           "perhaps you should check the basement for pods :)")))

jointData$weight <- runif(nrow(jointData), min=.98,max=1.02)
weightedFits <- c()

tdir <- paste0(tempdir(), "/")

write.cbor(mxData(jointData, "raw"), file=paste0(tdir, "joint-01.cbor"))
write.cbor(mxData(jointData[1:200,], "raw"), file=paste0(tdir, "joint-02.cbor"))
write.cbor(mxData(jointData, "raw"), file=paste0(tdir, "joint-03.cbor"))
write.cbor(mxData(jointData[1:200,], "raw"), file=paste0(tdir, "joint-04.cbor"))

jointModel1 <- mxModel("ContinuousOrdinalData",
	mxData(jointData[1:200,], "raw"),
	loadings, resid, means, thresh,
	mxAlgebra(t(L) %*% L + U, name="C"),
	mxFitFunctionML(),
	mxExpectationNormal("C", "M",
		dimnames=names(jointData)[1:5],
		thresholds="T",
		threshnames=c("z2", "z4", "z5")),
	mxComputeLoop(list(
		mxComputeLoadData("ContinuousOrdinalData", paste0(tdir, "joint-%02d.cbor")),
		mxComputeGradientDescent(),
		mxComputeCheckpoint("C", path="loadData.csv", toReturn=TRUE)
	), i=1:4))

jointModel1 <- mxRun(jointModel1)

log <- jointModel1$compute$steps[[3]]$log

dlog <- read.table("loadData.csv", sep="\t",
                   header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

omxCheckEquals(colnames(dlog), colnames(log))
omxCheckEquals(nrow(dlog), nrow(log))

for (cx in 1:ncol(log)) {
	cn <- colnames(dlog)[cx]
  if (cn == 'timestamp' || cn == 'statusCode') next
  omxCheckCloseEnough(log[[cx]], dlog[[cx]], 1e-9)
}

omxCheckCloseEnough(coef(jointModel1),
                    c(log[nrow(log), names(coef(jointModel1))]), 1e-9)
omxCheckCloseEnough(c(jointModel1$C$result),
                    c(log[nrow(log), grepl("ContinuousOrdinalData.C", colnames(log), fixed=TRUE)]), 1e-9)

ignCols <- -match(c('MxComputeLoop1','timestamp','iterations',
        'OpenMxEvals', 'statusCode'),
      colnames(log))
omxCheckCloseEnough(log[1,ignCols], log[3,ignCols], 1e-3)
omxCheckCloseEnough(log[2,ignCols], log[4,ignCols], 1e-3)
omxCheckTrue(max(abs(log[1,ignCols]- log[2,ignCols])) > 1)
