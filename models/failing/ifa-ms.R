library(OpenMx)
library(rpf)

m2.data <- suppressWarnings(try(read.table("models/failing/data/ms-data.csv"), silent=TRUE))
if (is(m2.data, "try-error")) m2.data <- read.table("data/ms-data.csv")
m2.data[m2.data==-9] <- NA
m2.data <- m2.data + 1

m2.data <- data.frame(lapply(m2.data, mxFactor, levels=1:5))

gpcm <- function(outcomes) {
  rpf.nrm(outcomes, T.c=lower.tri(diag(outcomes-1),TRUE) * -1)
  #   rpf.nrm(outcomes, T.c=diag(outcomes-1))
}

m2.spec <- list()
m2.spec[1:22] <- gpcm(5)
m2.spec[2] <- gpcm(4)
m2.spec[5] <- gpcm(3)
m2.spec[6] <- gpcm(4)
m2.spec[13:14] <- gpcm(4)

m2.numItems <- length(m2.spec)
m2.maxParam <-max(sapply(m2.spec, rpf.numParam))

ip.mat <- mxMatrix(name="ItemParam", nrow=m2.maxParam, ncol=m2.numItems,
                   values=c(1, 1, rep(0, m2.maxParam-2)), free=FALSE)
ip.mat@labels[1,] <- 'a1'
ip.mat@free[1,] <- TRUE
for (ix in 1:m2.numItems) {
  thr <- m2.spec[[ix]]@outcomes - 1
  ip.mat@free[(2+thr):(1+2*thr), ix] <- TRUE
}

#  m2.fmfit <- read.flexmirt("~/2012/sy/fm/ms-prm.txt")
#  ip.mat@values <- m2.fmfit$G1$param

eip.mat <- mxAlgebra(ItemParam, name="EItemParam")

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

m2 <- mxModel(model="m2", eip.mat, m.mat, cov.mat, ip.mat,
              mxData(observed=m2.data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec=m2.spec,
                                ItemParam="ItemParam", EItemParam="EItemParam"),
              mxFitFunctionML(),
              mxComputeIterate(steps=list(
                mxComputeOnce('expectation', context='EM'),
#                mxComputeGradientDescent(free.set='ItemParam', useGradient=TRUE, verbose=1L),
                mxComputeNewtonRaphson(free.set='ItemParam'),
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction', adjustStart=TRUE, free.set=c("mean", "cov"))
              )))
#  m2 <- mxOption(m2, "Number of Threads", 1)
m2 <- mxRun(m2)
omxCheckCloseEnough(m2@fitfunction@result, 50510, 5)  # unstable
