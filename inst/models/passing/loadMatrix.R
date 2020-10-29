library(OpenMx)
library(mvtnorm)

Cov <- rWishart(4, 20, toeplitz(c(2,1)/20))
write.table(t(apply(Cov, 3, vech)),
            file="cov.csv", col.names=FALSE, row.names=FALSE)
Mean <- rmvnorm(4, sigma = diag(2))
write.table(Mean, file="mean.csv", col.names=FALSE, row.names=FALSE)

m1 <- mxModel(
  "test1",
  mxMatrix("Full", 1,2, values=0,       name="mean"),
  mxMatrix("Symm", 2,2, values=diag(2), name="cov"),
  mxMatrix("Full", 1,2, values=-1,      name="lbound"),
  mxMatrix("Full", 1,2, values=1,       name="ubound"),
  mxAlgebra(omxMnor(cov,mean,lbound,ubound), name="area"),
  mxFitFunctionAlgebra("area"),
  mxComputeLoop(list(
    mxComputeLoadMatrix(c('mean', 'cov'),
	    path=c('mean.csv', 'cov.csv')),
    mxComputeOnce('fitfunction', 'fit'),
    mxComputeCheckpoint(path="loadMatrix.csv")
  ), i=1:4))

m1 <- mxRun(m1)

log <- read.table("loadMatrix.csv", sep="\t",
                   header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

e1 <- sapply(1:4, function (x) omxMnor(covariance = Cov[,,x], means = Mean[x,],
                                 lbound = c(-1,-1), ubound = c(1,1)))

omxCheckCloseEnough(log$objective, e1, 1e-4)

#---------

df <- as.data.frame(cbind(Mean, t(apply(Cov, 3, vech))))

plan <- mxComputeLoop(list(
  mxComputeLoadMatrix(c('mean', 'cov'), observed=df, method = 'data.frame'),
  mxComputeOnce('fitfunction', 'fit'),
  CP=mxComputeCheckpoint(toReturn=TRUE)), maxIter = 4L)

m1 <- mxRun(mxModel(m1, plan))

omxCheckCloseEnough(e1, m1$compute$steps$CP$log$objective, 1e-4)

#---------

m2 <- mxModel(m1,   mxComputeLoop(list(
    mxComputeLoadMatrix(c('mean', 'cov'),
	    path=c('mean.csv', 'cov.csv')),
    mxComputeOnce('fitfunction', 'fit'),
    mxComputeCheckpoint(path="loadMatrix.csv")
), i=3:4))

m2 <- mxRun(m2)

log <- read.table("loadMatrix.csv", sep="\t",
                   header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

omxCheckCloseEnough(log$objective, e1[3:4], 1e-4)

m3 <- mxModel(m1,   mxComputeLoop(list(
    mxComputeLoadMatrix(c('mean', 'cov'),
	    path=c('mean.csv', 'cov.csv')),
    mxComputeOnce('fitfunction', 'fit'),
    mxComputeCheckpoint(path="loadMatrix.csv")
), i=2:1))

omxCheckError(mxRun(m3),
	"MxComputeLoadMatrix: at line 3, cannot seek backwards to line 1")
