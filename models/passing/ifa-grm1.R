#options(error = browser)
library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 14
spec <- vector("list", numItems)
for (ix in 1:numItems) { spec[[ix]] <- rpf.grm(outcomes=sample(2:7, 1)) }
correct <- lapply(spec, rpf.rparam, version=1)

ability <- rnorm(500)
data <- rpf.sample(ability, spec, correct)

ip.mat <- mxMatrix(name="itemParam", nrow=max(sapply(correct, length)), ncol=numItems)
colnames(ip.mat) <- colnames(data)
rownames(ip.mat) <- c('f1', paste('b', 1:(nrow(ip.mat)-1), sep=""))
correct.mat <- ip.mat$values
for (ix in 1:numItems) {
  len <- length(correct[[ix]])
  ip.mat$free[1:len,ix] <- TRUE
  ip.mat$values[1:len,ix] <- rpf.rparam(spec[[ix]], version=1)
  correct.mat[1:len,ix] <- correct[[ix]]
}
ip.mat$values[!ip.mat$free] <- NA
correct.mat[!ip.mat$free] <- NA

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
rownames(m.mat) <- "f1"
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)
dimnames(cov.mat) <- list("f1", "f1")

m2 <- mxModel(model="grm1", ip.mat, m.mat, cov.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=spec, ItemParam="itemParam",
                mean="mean", cov="cov", qpoints=31, debugInternal=TRUE),
              mxFitFunctionML(),
	      mxComputeOnce('expectation', 'scores'))
middle <- mxRun(m2)
omxCheckCloseEnough(sum(middle$expectation$debug$patternLikelihood), -9742.31, .1)
omxCheckCloseEnough(fivenum(middle$expectation$debug$patternLikelihood),
                    c(-34.98313, -22.50933, -19.59691, -16.79153, -6.51683), .001)
omxCheckCloseEnough(sum(middle$expectation$debug$em.expected), 7000, .01)
omxCheckCloseEnough(fivenum(middle$expectation$debug$em.expected),
                    c(0, 0, 0.00451, 1.61901, 80.99209), .01)

testDeriv <- mxModel(m2,
	      mxComputeSequence(list(
		  mxComputeOnce('expectation', 'scores'),
		  mxComputeOnce('fitfunction', c('fit', 'gradient', 'hessian', 'ihessian')),
      mxComputeReportDeriv()
				)))
testDeriv <- mxRun(testDeriv)
omxCheckCloseEnough(testDeriv$fitfunction$result, 2*9399.954, .01)
omxCheckCloseEnough(fivenum(testDeriv$output$gradient),
                    2*c(-14424.48407, -62.52714, -2.51876, 71.87544, 14651.19963), .01)
omxCheckCloseEnough(fivenum(testDeriv$output$hessian[testDeriv$output$hessian != 0]),
                    2*c(-1038404.94356, -20.82415, 0.06734, 53.01333, 1038503.00369 ), .01)
omxCheckCloseEnough(max(abs(solve(testDeriv$output$hessian) - testDeriv$output$ihessian)), 0, .001)

plan <- mxComputeSequence(list(mxComputeEM('expectation', 'scores',
                                           mxComputeNewtonRaphson()),
                               mxComputeOnce('fitfunction', 'information', "meat"),
                               mxComputeHessianQuality(),
                               mxComputeStandardError()))

m2 <- mxModel(m2,
              mxExpectationBA81(
                ItemSpec=spec, ItemParam="itemParam",
                mean="mean", cov="cov",
                qpoints=31,
                scores="full"),
	      plan)
				  
m2 <- mxRun(m2)

emstat <- m2$compute$steps[[1]]$output
omxCheckCloseEnough(emstat$EMcycles, 14, 1)
omxCheckCloseEnough(emstat$totalMstep, 51, 5)

#print(m2$matrices$itemParam$values)
#print(correct.mat)
omxCheckCloseEnough(m2$fitfunction$result, 13969.747, .01)
got <- cor(c(m2$matrices$itemParam$values[ip.mat$free]),
           c(correct.mat[ip.mat$free]))
omxCheckCloseEnough(got, .993, .01)
scores <- m2$expectation$output$scores
omxCheckCloseEnough(scores[1:5,1], c(0.81609, 0.74994, -0.83515, 0.79766, 0.16879), .001)
omxCheckCloseEnough(scores[1:5,2], c(0.43522, 0.44211, 0.4686, 0.43515, 0.3842), .001)
omxCheckCloseEnough(sum(abs(scores[,1] - ability) < 1*scores[,2])/500, .714, .01)
omxCheckCloseEnough(sum(abs(scores[,1] - ability) < 2*scores[,2])/500, .95, .01)

omxCheckTrue(m2$output$infoDefinite)
omxCheckCloseEnough(c(m2$output$conditionNumber), 658.58, 1)

#cat(deparse(round(m2$output$standardErrors,3)))

se <- c(0.152, 0.114, 0.115, 0.261, 0.163, 0.121, 0.142,  0.109, 0.113, 0.098,
        0.185, 0.168, 0.161, 0.132, 0.139, 0.154,  0.126, 0.113, 0.215, 0.107,
        0.111, 0.149, 0.124, 0.12, 0.127,  0.186, 0.135, 0.134, 0.134, 0.145,
        0.147, 0.125, 0.15, 0.12,  0.12, 0.122, 0.123, 0.152, 0.164, 0.129, 0.117,
        0.104, 0.109,  0.139, 0.103, 0.104, 0.106, 0.109, 0.116, 0.135, 0.163, 0.115 )
  
omxCheckCloseEnough(c(m2$output$standardErrors), se, .01)

i2 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction', 'information', "sandwich"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i2 <- mxRun(i2, silent=TRUE)

omxCheckCloseEnough(i2$output$conditionNumber, 662, 1)
#cat(deparse(round(i2$output$standardErrors,3)))
swse <- c(0.143, 0.11, 0.11, 0.238, 0.149, 0.125, 0.134, 0.106,  0.108, 0.094,
          0.169, 0.165, 0.154, 0.128, 0.13, 0.148, 0.121,  0.119, 0.199, 0.102,
          0.106, 0.144, 0.117, 0.115, 0.122, 0.177,  0.132, 0.13, 0.129, 0.136,
          0.141, 0.134, 0.145, 0.116, 0.112,  0.116, 0.116, 0.151, 0.162, 0.124,
          0.115, 0.097, 0.104, 0.125,  0.098, 0.099, 0.104, 0.107, 0.111, 0.139, 0.156, 0.11)
omxCheckCloseEnough(c(i2$output$standardErrors), swse, .001)
