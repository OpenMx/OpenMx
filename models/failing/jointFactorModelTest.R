library(OpenMx)
library(mvtnorm)
jointDataCov <- rbind(  c( 2, .3, .6,  .9, .6),
                        c(.3,  4, .4,  .3, .9),
                        c(.6, .4,  5,  1.2, .0),
                        c(.9, .3, 1.2,  1,  .4),
                        c(.6, .9, .0,  .4,  3))
                        
jointData <- data.frame(rmvnorm(2, c(8, 0, 2, 0, 7, 0), jointDataCov))
names(jointData) <- paste("z", 1:6, sep="")
jointData[,2] <- as.integer(jointData[,2] > 0)
tester = jointData[,4]
jointData[, 4] <- 0
jointData[(tester > -.35), 4] <- 1
jointData[(tester > .15), 4] <- 2
jointData[(tester > .75), 4] <- 3

tester = jointData[,5]
jointData[,5] <- 0
jointData[(tester > -.8) ,5] <- 1
jointData[(tester > -.3),5] <- 2

# specify ordinal columns as ordered factors
jointData[,c(2,4,5)] <- mxFactor(jointData[,c(2,4,5)], 
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))
	
satCov <- mxMatrix("Symm", 5, 5,
	free=TRUE, values=diag(5), name="C")
satCov@free[2,2] <- FALSE
satCov@free[4,4] <- FALSE
satCov@free[5,5] <- FALSE

loadings <- mxMatrix("Full", 1, 5,
	free=TRUE, values=1, name="L")
	
resid <- mxMatrix("Diag", 5, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=.5, name="U")
	
means <- mxMatrix("Full", 1, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=0, name="M")
	
thresh <- mxMatrix("Full", 3, 3, FALSE, 0, name="T")

thresh@free[,1] <- c(TRUE, FALSE, FALSE)
thresh@values[,1] <- c(0, NA, NA)
thresh@labels[,1] <- c("z2t1", NA, NA)

thresh@free[,2] <- TRUE
thresh@values[,2] <- c(-1, 0, 1)
thresh@labels[,2] <- c("z4t1", "z4t2", "z4t3")

thresh@free[,3] <- c(TRUE, TRUE, FALSE)
thresh@values[,3] <- c(-1, 1, NA)
thresh@labels[,3] <- c("z5t1", "z5t2", NA)
	
# 
jointModel1 <- mxModel("ContinuousOrdinalData",
	mxData(jointData, "raw"),
	loadings, resid, means, thresh,
	mxAlgebra(t(L) %*% L + U, name="C"),
	mxFIMLObjective("C", "M",
		dimnames=names(jointData),
		thresholds="T",
		threshnames=c("z2", "z4", "z5"))
	)
	
jointResults1 <- mxRun(jointModel1, unsafe=TRUE)

jointModel2 <- mxModel("ContinuousOrdinalData",
  mxData(jointData, "raw"),
  satCov, means, thresh,
  mxFIMLObjective("C", "M",
      dimnames=names(jointData),
      thresholds="T",
      threshnames=c("z2", "z4", "z5"))
  )
  
jointResults2 <- mxRun(jointModel2, unsafe=TRUE)