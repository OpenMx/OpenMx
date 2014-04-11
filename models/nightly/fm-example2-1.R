# flexMIRT 1.88/2.0 Example 2-1 "2PL Calibration Syntax"

library(OpenMx)
library(rpf)

Scale <- -2
mxOption(NULL, 'loglikelihoodScale', Scale)

set.seed(1)
m2.data <- suppressWarnings(try(read.table("models/nightly/data/g341-19.dat"), silent=TRUE))
if (is(m2.data, "try-error")) m2.data <- read.table("data/g341-19.dat")
m2.data <- m2.data + 1

m2.spec <- list()
m2.spec[1:12] <- rpf.grm(outcomes=2)
m2.numItems <- length(m2.spec)

for (c in 1:m2.numItems) {
  m2.data[[c]] <- mxFactor(m2.data[[c]], levels=1:m2.spec[[c]]$outcomes)
}

m2.maxParam <-max(sapply(m2.spec, rpf.numParam))

ip.mat <- mxMatrix(name="ItemParam", nrow=m2.maxParam, ncol=m2.numItems,
                   values=c(1, 0), free=TRUE)

#  m2.fmfit <- read.flexmirt("~/2012/sy/fm/ms-rasch-prm.txt")
# cat(deparse(round(m2.fmfit$G1$param,6)))
#  ip.mat$values <- m2.fmfit$G1$param

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

m2 <- mxModel(model="m2", m.mat, cov.mat, ip.mat,
              mxData(observed=m2.data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec=m2.spec,
                                ItemParam="ItemParam"),
              mxFitFunctionML(),
	      mxComputeSequence(steps=list(
				    mxComputeEM('expectation', 'scores',
				                mxComputeNewtonRaphson(freeSet='ItemParam'),
				                mxComputeOnce('fitfunction', 'fit'),
				                information=TRUE, semDebug=TRUE, info.method="hessian",
						infoArgs=list(fitfunction=c('fitfunction'))),
				    mxComputeStandardError(),
				    mxComputeHessianQuality(),
            mxComputeReportDeriv())))

#  m2 <- mxOption(m2, "Number of Threads", 1)
m2 <- mxRun(m2, silent=TRUE)

emstat <- m2$compute$steps[[1]]$output
omxCheckCloseEnough(emstat$EMcycles, 17, 1)
omxCheckCloseEnough(emstat$totalMstep, 53, 10)
omxCheckCloseEnough(emstat$semProbeCount, 96, 3)

omxCheckCloseEnough(m2$output$minimum, Scale * 33408.05/-2, .01)
omxCheckTrue(m2$output$infoDefinite)
omxCheckCloseEnough(m2$output$conditionNumber, 51, 1)
#cat(deparse(round(c(m2$output$standardErrors),3)))

se <- c(0.071, 0.047, 0.109, 0.113, 0.062, 0.043, 0.073, 0.052, 0.063,
        0.044, 0.098, 0.086, 0.181, 0.233, 0.075, 0.061, 0.153, 0.166,
        0.102, 0.096, 0.064, 0.044, 0.069, 0.045)
omxCheckCloseEnough(c(m2$output$standardErrors), se, .01)
#max(abs(c(m2$output$standardErrors) - se))

m3 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation', 'scores'),
                mxComputeOnce('fitfunction', 'information', "hessian"),
                mxComputeReportDeriv())))
m3 <- mxRun(m3, silent=TRUE)

m5 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction', 'fit'),
                mxComputeNumericDeriv(parallel=FALSE, iterations=2L),
                mxComputeStandardError(),
                mxComputeHessianQuality(),
                mxComputeReportDeriv())))
m5 <- mxRun(m5, silent=TRUE)
omxCheckTrue(m5$output$infoDefinite)
omxCheckCloseEnough(m5$output$conditionNumber, 51, 1)

if (0) {
  probe <- function(pt) {
    ip.mat$values[,] <- pt
    m6 <-mxModel(m2, ip.mat,
                 mxComputeSequence(steps=list(
                   mxComputeOnce('expectation'),
                   mxComputeOnce('fitfunction', 'fit'))))
    m6 <- mxRun(m6, silent=TRUE)
    fit <- m6$output$minimum
    fit
  }
  require(numDeriv)
  H.nd <- hessian(probe, m2$matrices$ItemParam$values, method.args=list(r=2))
}

quantifyAsymmetry <- function(info) {
  sym1 <- (info + t(info))/2
  sym2 <- try(chol(solve(sym1)), silent=TRUE)
  if (inherits(sym2, "try-error")) return(NA)
  asymV <- (info - t(info))/2
  max(svd(sym2 %*% asymV %*% sym2, 0, 0)$d)
}
  
if(1) {
  iinfo <- m2$compute$steps[[1]]$debug$inputInfo
#  print(iinfo[1:5,1:5])
  dm  <- m2$compute$steps[[1]]$debug$rateMatrix
#  print(dm[1:5,1:5])
  #  dm1 <- (dm + t(dm)) / 2
  dm1 <- t(dm)   # without averaging with transpose
  dm2 <- diag(dim(dm)[1]) - dm1
  H <- m3$output$hessian
  tmp <- (dm2 %*% H)
  emHess <- solve(tmp)
  tmp2 <- (tmp + t(tmp))/2
  emHess2 <- solve(tmp2)
  kappa(emHess2, exact=TRUE)
  omxCheckCloseEnough(max(abs(m2$output$ihessian - emHess2)), 0, 1e-4)
  #  sv <- svd(emHess)$d
#  max(sv)/min(sv)
  omxCheckCloseEnough(kappa(emHess, exact=TRUE), 51, 1)
#  kappa(m2$output$ihessian, exact=TRUE)
  
  omxCheckCloseEnough(max(abs(diag(emHess) - diag(solve(m5$output$hessian)))), 0, .001)
  omxCheckCloseEnough(quantifyAsymmetry(emHess), .071, .05)
  omxCheckCloseEnough(quantifyAsymmetry(emHess2), 0, 1e-6)
  #hist(abs(diag(emHess) - diag(solve(m5$output$hessian))))
  
  omxCheckCloseEnough(max(sqrt(abs(Scale)*diag(emHess)) - c(m2$output$standardErrors)), 0, 5e-5)
  #print(m2$matrices$ItemParam$values - fmfit)
}

print(m2$output$backendTime)
