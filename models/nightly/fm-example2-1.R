# flexMIRT 1.88/2.0 Example 2-1 "2PL Calibration Syntax"

library(OpenMx)
library(rpf)

set.seed(1)
m2.data <- suppressWarnings(try(read.table("models/nightly/data/g341-19.dat"), silent=TRUE))
if (is(m2.data, "try-error")) m2.data <- read.table("data/g341-19.dat")
m2.data <- m2.data + 1

m2.spec <- list()
m2.spec[1:12] <- rpf.grm(outcomes=2)
m2.numItems <- length(m2.spec)

for (c in 1:m2.numItems) {
  m2.data[[c]] <- mxFactor(m2.data[[c]], levels=1:m2.spec[[c]]@outcomes)
}

m2.maxParam <-max(sapply(m2.spec, rpf.numParam))

ip.mat <- mxMatrix(name="ItemParam", nrow=m2.maxParam, ncol=m2.numItems,
                   values=c(1, 0), free=TRUE)

#  m2.fmfit <- read.flexmirt("~/2012/sy/fm/ms-rasch-prm.txt")
# cat(deparse(round(m2.fmfit$G1$param,6)))
#  ip.mat@values <- m2.fmfit$G1$param

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

m2 <- mxModel(model="m2", m.mat, cov.mat, ip.mat,
              mxData(observed=m2.data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec=m2.spec,
                                ItemParam="ItemParam"),
              mxFitFunctionML(),
              mxComputeEM('expectation',
			  mxComputeNewtonRaphson(free.set='ItemParam'),
			  mxComputeOnce('fitfunction', fit=TRUE, free.set=c("mean", "cov"))))

#  m2 <- mxOption(m2, "Number of Threads", 1)
m2 <- mxRun(m2, silent=TRUE)
omxCheckCloseEnough(m2@fitfunction@result, 33408.05, .01)

#print(m2@matrices$ItemParam@values - fmfit)
print(m2@output$backendTime)
