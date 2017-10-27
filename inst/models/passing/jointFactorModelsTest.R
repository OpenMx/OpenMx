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

for (strat in c('auto', 'ordinal', 'continuous')) {
	print(paste('***',strat,'***'))
					# run factor and saturated models
	jointModel1 <- mxModel("ContinuousOrdinalData",
			       mxData(jointData, "raw"),
			       loadings, resid, means, thresh,
			       mxAlgebra(t(L) %*% L + U, name="C"),
			       mxFitFunctionML(jointConditionOn=strat),
			       mxExpectationNormal("C", "M",
						   dimnames=names(jointData)[1:5],
						   thresholds="T",
						   threshnames=c("z2", "z4", "z5"))
			       )
	
	jointResults1 <- mxRun(jointModel1, suppressWarnings = TRUE)
	#summary(jointResults1)
	print(jointResults1$fitfunction$info)

	#cat(deparse(round(coef(jointResults1),3)))
	c1 <- c(0.609, 0.579, 0.657, 0.606, 0.165, 0.551, 0.499,
		7.978, 2.069, 0.059, -0.387, 0.116, 0.815, -0.633, -0.285)
	print(max(abs(coef(jointResults1) - c1)))
	omxCheckCloseEnough(coef(jointResults1), c1, 0.005)

	omxCheckCloseEnough(jointResults1$output$Minus2LogLikelihood, 2683.071, 0.2)

	jointModel2 <- mxModel("ContinuousOrdinalData",
			       mxData(jointData, "raw"),
			       satCov, means, thresh,
			       mxFitFunctionML(jointConditionOn=strat),
			       mxExpectationNormal("C", "M",
						   dimnames=names(jointData)[1:5],
						   thresholds="T",
						   threshnames=c("z2", "z4", "z5"))
			       )
	
	jointResults2 <- mxRun(jointModel2, suppressWarnings = TRUE)
	#summary(jointResults2)

	#cat(deparse(round(coef(jointResults2),3)))
	c2 <- c(0.923, 0.374, 0.436, 0.367, 0.933, 0.354, 0.48, 0.43,  0.139,
		0.266, 0.083, 0.201, 7.978, 2.07, 0.065, -0.413, 0.127,  0.874, -0.872, -0.394)
	print(max(abs(coef(jointResults2) - c2)))
	omxCheckCloseEnough(coef(jointResults2), c2, 0.005)
	
	omxCheckCloseEnough(jointResults2$output$Minus2LogLikelihood, 2674.235, 0.1)

	jointModel3 <- mxModel(jointResults2)
	jointModel3$data$weight <- 'weight'
	weightedFit <- mxRun(mxModel(jointModel3, mxComputeOnce('fitfunction','fit')))
	weightedFits <- c(weightedFits, weightedFit$output$fit)
}

omxCheckCloseEnough(diff(weightedFits), rep(0,2), .04)
