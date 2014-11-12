require(OpenMx)

# get data
jointData <- suppressWarnings(try(read.table("models/passing/data/jointdata.txt", header=TRUE), silent=TRUE))
jointData <- read.table("data/jointdata.txt", header=TRUE)

# specify ordinal columns as ordered factors
jointData[,c(2,4,5)] <- mxFactor(jointData[,c(2,4,5)], 
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))
	
satCov <- mxMatrix("Symm", 5, 5,
	free=TRUE, values=diag(5), name="C")
satCov$free[2,2] <- FALSE
satCov$free[4,4] <- FALSE
satCov$free[5,5] <- FALSE

loadings <- mxMatrix("Full", 1, 5,
	free=TRUE, values=1, name="L")
	
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

# run factor and saturated models
jointModel1 <- mxModel("ContinuousOrdinalData",
	mxData(jointData, "raw"),
	loadings, resid, means, thresh,
	mxAlgebra(t(L) %*% L + U, name="C"),
	mxFitFunctionML(),mxExpectationNormal("C", "M",
		dimnames=names(jointData),
		thresholds="T",
		threshnames=c("z2", "z4", "z5"))
	)
	
jointResults1 <- mxRun(jointModel1, suppressWarnings=TRUE)

jointModel2 <- mxModel("ContinuousOrdinalData",
  mxData(jointData, "raw"),
  satCov, means, thresh,
  mxFitFunctionML(),mxExpectationNormal("C", "M",
      dimnames=names(jointData),
      thresholds="T",
      threshnames=c("z2", "z4", "z5"))
  )
  
jointResults2 <- mxRun(jointModel2, suppressWarnings=TRUE)

# check that the likelihoods are the same as originally reported
omxCheckCloseEnough(jointResults1$output$Minus2LogLikelihood, 2683.071, 0.01)
omxCheckCloseEnough(jointResults2$output$Minus2LogLikelihood, 2674.235, 0.01)
