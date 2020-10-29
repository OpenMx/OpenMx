library(OpenMx)

if (mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")
mxOption(NULL, "major iterations", 3000)

data("jointdata", package ="OpenMx", verbose= TRUE)

jointdata[,c("z2","z4","z5")] <- mxFactor(jointdata[,c("z2","z4","z5")], 
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

satCov <- mxMatrix(name="C", "Symm", 5, 5, free=TRUE, values=diag(5))
satCov$free[2,2] <- FALSE
satCov$free[4,4] <- FALSE
satCov$free[5,5] <- FALSE

loadings <- mxMatrix(name="L", "Full", 1, 5, free=TRUE, values=1, lbound=0)
loadings$ubound[1,] <- 2
	
resid <- mxMatrix(name="U", "Diag", 5, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=.5)
	
means <- mxMatrix(name="M", "Full", 1, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=0)
	
thresh <- mxMatrix(name="T", "Full", 3, 3, FALSE, 0)

thresh$free[,1]   <- c(TRUE, FALSE, FALSE)
thresh$values[,1] <- c(0, NA, NA)
thresh$labels[,1] <- c("z2t1", NA, NA)

thresh$free[,2]   <- TRUE
thresh$values[,2] <- c(-1, 0, 1)
thresh$labels[,2] <- c("z4t1", "z4t2", "z4t3")

thresh$free[,3]   <- c(TRUE, TRUE, FALSE)
thresh$values[,3] <- c(-1, 1, NA)
thresh$labels[,3] <- c("z5t1", "z5t2", NA)
colnames(thresh)  <- paste0('z', c(2,4,5))

jm1 <- mxModel("ContinuousOrdinalData", loadings, resid, means, thresh,
			mxAlgebra(t(L) %*% L + U, name="C"),
			mxFitFunctionWLS(),
			mxExpectationNormal("C", "M",
				dimnames=names(jointdata),
				thresholds="T",
				threshnames=c("z2", "z4", "z5")),
			mxDataWLS(jointdata, "WLS"),
			mxCI('L', interval = .8)
)

jm1 <- mxRun(jm1, intervals = TRUE)

jm2 <- mxModel(jm1, mxData(jointdata, 'raw'), mxFitFunctionML())
jm2 <- mxRun(jm2, intervals=TRUE)

ci1 <- jm1$output$confidenceIntervals
ci2 <- jm2$output$confidenceIntervals
goodEntry <- !is.na(ci1) & !is.na(ci2)

omxCheckCloseEnough(sum(goodEntry), 15, 2)

print(max(abs(ci1[goodEntry] - ci2[goodEntry])))
omxCheckCloseEnough(ci1[goodEntry], ci2[goodEntry], .09)
omxCheckCloseEnough(median(abs(ci1[goodEntry] - ci2[goodEntry])), 0, .02)

# Just FYI... WLS is 70 times quicker even with only n= 175 
# and just 2,3, and 4 levels in the ordinal vars!
# umx_time(jm1) # ContinuousOrdinalData: 00.12 seconds
# umx_time(jm2) # ContinuousOrdinalData: 08.18 seconds
