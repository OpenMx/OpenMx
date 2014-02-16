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
				                information=TRUE, semDebug=TRUE),  #, verbose=2L, semDebug=TRUE
				    mxComputeStandardError(),
				    mxComputeHessianQuality())))

#  m2 <- mxOption(m2, "Number of Threads", 1)
m2 <- mxRun(m2, silent=TRUE)
omxCheckCloseEnough(m2@output$minimum, 33408.05, .01)
#omxCheckTrue(m2@output$infoDefinite)
#omxCheckCloseEnough(m2@output$conditionNumber, 61, 1)

#cat(deparse(round(c(m2@output$standardErrors),3)))
# se <- c(0.145, 0.091, 0.215, 0.174, 0.118, 0.077, 0.145, 0.093,
#         0.118,  0.077, 0.208, 0.143, 0.408, 0.353, 0.144, 0.097,
#         0.369, 0.269,  0.205, 0.148, 0.121, 0.07, 0.137, 0.073)
# omxCheckCloseEnough(c(m2@output$standardErrors), se, .04)
#max(abs(c(m2@output$standardErrors) - se))

m3 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation', context="EM"),
                mxComputeOnce('fitfunction', information=TRUE, info.method="meat"))))
m3 <- mxRun(m3, silent=TRUE)

m4 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction', information=TRUE, info.method="meat"))))
m4 <- mxRun(m4, silent=TRUE)
omxCheckCloseEnough(max(abs(m3@output$hessian - m4@output$hessian)), 0, 1e-12)

m5 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation', context=""),
                mxComputeNumericDeriv(parallel=FALSE, iterations=2L),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
m5 <- mxRun(m5, silent=TRUE)
omxCheckTrue(m5@output$infoDefinite)
omxCheckCloseEnough(m5@output$conditionNumber, 51, 1)

if(0) {
  is.symm <- function(a) max(abs(a- t(a)))
  dm  <- m2@compute@steps[[1]]@debug$rateMatrix
  dm1 <- (dm + t(dm)) / 2
  is.symm(dm1)
  dm2 <- diag(dim(dm)[1]) - dm1
  is.symm(dm2)
  H <- m3@output$hessian
  is.symm(H)
  iH <- solve(H)
  is.symm(iH)
  is.symm(solve(dm2))
  emHess <- iH %*% solve(dm2)
  emHess1 <- (emHess + t(emHess))/2
  emHess2 <- solve(dm2 %*% H)
  is.symm(emHess2)
  eigen(emHess)$val
  kappa(emHess, exact=TRUE)
  kappa(emHess1, exact=TRUE)
  sqrt(2*diag(emHess)) - m2@output$standardErrors
  
  max(abs(diag(emHess) - diag(solve(m5@output$hessian))))
  
  #print(m2@matrices$ItemParam@values - fmfit)
  print(m2@output$backendTime)
}
