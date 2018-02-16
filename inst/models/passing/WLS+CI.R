library(OpenMx)

if (mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")
mxOption(NULL, "major iterations", 3000)

jointData <- suppressWarnings(try(read.table("models/passing/data/jointdata.txt", header=TRUE), silent=TRUE))
jointData <- read.table("data/jointdata.txt", header=TRUE)

jointData[,c(2,4,5)] <- mxFactor(jointData[,c(2,4,5)], 
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

satCov <- mxMatrix("Symm", 5, 5,
	free=TRUE, values=diag(5), name="C")
satCov$free[2,2] <- FALSE
satCov$free[4,4] <- FALSE
satCov$free[5,5] <- FALSE

loadings <- mxMatrix("Full", 1, 5,
	free=TRUE, values=1, name="L", lbound=0)
loadings$ubound[1,] <- 2
	
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
colnames(thresh) <- paste0('z', c(2,4,5))

jm1 <- mxModel("ContinuousOrdinalData",
				mxDataWLS(jointData, "WLS"),
				loadings, resid, means, thresh,
			mxAlgebra(t(L) %*% L + U, name="C"),
			mxFitFunctionWLS(),
			mxExpectationNormal("C", "M",
				dimnames=names(jointData),
				thresholds="T",
				threshnames=c("z2", "z4", "z5")),
			mxCI('L', interval=.8)
			)

jm1 <- mxRun(jm1, intervals = TRUE)

jm2 <- mxModel(jm1, mxData(jointData, 'raw'),
               mxFitFunctionML())
jm2 <- mxRun(jm2, intervals=TRUE)

ci1 <- jm1$output$confidenceIntervals
ci2 <- jm2$output$confidenceIntervals
goodEntry <- !is.na(ci1) & !is.na(ci2)

omxCheckCloseEnough(sum(goodEntry), 15, 2)

print(max(abs(ci1[goodEntry] - ci2[goodEntry])))
omxCheckCloseEnough(ci1[goodEntry], ci2[goodEntry], .08)
