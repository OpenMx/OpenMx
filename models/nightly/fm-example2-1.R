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
	      mxComputeSequence(steps=list(
				    mxComputeEM('expectation',
				                mxComputeNewtonRaphson(free.set='ItemParam'),
				                mxComputeOnce('fitfunction', fit=TRUE, free.set=c("mean", "cov")),
				                information=TRUE),  #, verbose=2L, semDebug=TRUE
				    mxComputeStandardError(),
				    mxComputeHessianQuality())))

#  m2 <- mxOption(m2, "Number of Threads", 1)
m2 <- mxRun(m2, silent=TRUE)
omxCheckCloseEnough(m2@output$minimum, 33408.05, .01)
omxCheckTrue(!m2@output$infoDefinite)  # something screwed up TODO

m3 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation', context=""),
                mxComputeEstimatedHessian(parallel=FALSE, iterations=2L),
                mxComputeHessianQuality())))
m3 <- mxRun(m3, silent=TRUE)
omxCheckTrue(m3@output$infoDefinite)
omxCheckCloseEnough(m3@output$conditionNumber, 51, 1)

#cat(deparse(round(c(m2@output$standardErrors),3)))
se <- c(0.101, 0.067, 0.168, 0.156, 0.088, 0.062, 0.104, 0.074, 0.088,
        0.063, 0.147, 0.118, 0.328, 0.336, 0.109, 0.085, 0.263, 0.235,
        0.156, 0.133, 0.09, 0.063, 0.098, 0.066)
omxCheckCloseEnough(c(m2@output$standardErrors), se, .01)

#print(m2@matrices$ItemParam@values - fmfit)
print(m2@output$backendTime)
