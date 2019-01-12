library(OpenMx)

set.seed(1)

jointData <- suppressWarnings(try(read.table("models/passing/data/jointdata.txt", header=TRUE), silent=TRUE))
jointData <- read.table("data/jointdata.txt", header=TRUE)
jointData[,c(2,4,5)] <- mxFactor(jointData[,c(2,4,5)], 
                                 levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

jointData$z1c <- with(jointData, z1 * .1 + rnorm(nrow(jointData)))

jointData$z2c <- with(jointData, rnorm(nrow(jointData), mean=unclass(z2)*.2))

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

# ------- LISREL

jl <- mxModel("JointLISREL", type="LISREL", thresh,
              manifestVars = list(endogenous=paste0('z',1:5)),
              latentVars = list(endogenous=c('G','z1c','z2c')),
              mxData(jointData, "raw", verbose=0L),
              mxPath('one', paste0('z', c(1,3)), free=TRUE, labels=c('z1','z3')),
              mxPath(paste0('z', c(1,3)), arrows=2, free=TRUE, values=.5),
              mxPath(paste0('z', c(2,4,5)), arrows=2, free=FALSE, values=.5),
              mxPath('G', arrows=2, values=1, free=FALSE),
              mxPath('G', paste0('z', 1:5), free=TRUE, values=1, lbound=0, ubound=4),
              mxPath('one', 'z1c', free=FALSE, labels="data.z1c"),
              mxPath('one', 'z2c', free=FALSE, labels="data.z2c"),
              mxPath('z1c', 'z1', labels="r1"),
              mxPath('z2c', 'z2', labels="r2"),
              mxFitFunctionWLS())

jl$expectation$thresholds <- 'T'
#jl$expectation$verbose <- 1L

jl <- mxRun(jl)

#----------- RAM

buildModel <- function(manifests=paste0('z', 1:5), type='WLS', slopes=TRUE) {
	jr <- mxModel("JointRAM", type="RAM", thresh,
		manifestVars = manifests,
		latentVars = c('G','z1c','z2c'),
		mxData(jointData, "raw", verbose=0L),
		mxPath('one', paste0('z', c(1,3)), free=TRUE, labels=c('z1','z3')),
		mxPath(paste0('z', c(1,3)), arrows=2, free=TRUE, values=.5),
		mxPath(paste0('z', c(2,4,5)), arrows=2, free=FALSE, values=.5),
		mxPath('G', arrows=2, values=1, free=FALSE),
		mxPath('G', paste0('z', 1:5), free=TRUE, values=1, lbound=0),
		mxFitFunctionWLS(type))
	if (slopes) {
		jr <- mxModel(jr,
			mxPath('one', 'z1c', free=FALSE, labels="data.z1c"),
			mxPath('one', 'z2c', free=FALSE, labels="data.z2c"),
			mxPath('z1c', 'z1', labels="r1"),
			mxPath('z2c', 'z2', labels="r2"))
	}
	jr$expectation$thresholds <- 'T'
	jr
}

jm1 <- buildModel()

jm1 <- mxRun(jm1)
summary(jm1)

print(max(abs(coef(jl) - coef(jm1))))
omxCheckCloseEnough(coef(jl), coef(jm1), 2e-5)

# tmp <- c(jm1$output$standardErrors)
# names(tmp) <- c()
# cat(deparse(round(tmp,4)))

param <-  c(0.5808, 0.5903, 0.6526, 0.6066, 0.1745, 0.0871, 0.0504, 0.5476,  0.4639, 7.8323,
            2.0759, 0.0785, -0.3664, 0.1271, 0.7919, -0.6475,  -0.296)
omxCheckCloseEnough(coef(jm1), param, 1e-3)

param.se <- c(0.0613, 0.1056, 0.0684, 0.0942, 0.0665, 0.0541, 0.0648, 0.0559,  0.0644, 0.1054,
              0.0593, 0.0819, 0.0777, 0.0726, 0.0922, 0.0655,  0.0585)
omxCheckCloseEnough(c(jm1$output$standardErrors), param.se, 1e-3)

jm2 <- mxModel(jm1, mxFitFunctionML())
jm2 <- mxRun(jm2)
summary(jm2)

estCmp <- cbind(coef(jm2), coef(jm1))
omxCheckCloseEnough(cor(estCmp)[2,1], 1, 1e-4)

seCmp <- cbind(jm2$output$standardErrors, jm1$output$standardErrors)
omxCheckCloseEnough(cor(seCmp)[2,1], 1, .18)

# ------- Test permutation code

for (slopes in c(TRUE,FALSE)) {
	for (type in c('WLS','DWLS','ULS')) {
		jm3 <- buildModel(type=type, slopes=slopes)
		jm3$data$verbose <- 1L
		jm3 <- mxRun(jm3)
		jm4 <- mxModel(buildModel(paste0('z', 5:1), type=type, slopes=slopes), jm3$data)
		jm4 <- mxRun(jm4)
	}
}
