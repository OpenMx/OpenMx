#options(error = utils::recover)
library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 14
spec <- vector("list", numItems)
for (ix in 1:numItems) { spec[[ix]] <- rpf.grm(outcomes=sample(2:7, 1)) }
correct <- lapply(spec, rpf.rparam)

ability <- rnorm(500)
data <- rpf.sample(ability, spec, correct)

ip.mat <- mxMatrix(name="itemParam", nrow=max(sapply(correct, length)), ncol=numItems)
correct.mat <- ip.mat@values
for (ix in 1:numItems) {
  len <- length(correct[[ix]])
  ip.mat@free[1:len,ix] <- TRUE
  ip.mat@values[1:len,ix] <- rpf.rparam(spec[[ix]])
  correct.mat[1:len,ix] <- correct[[ix]]
}
ip.mat@values[!ip.mat@free] <- NA
correct.mat[!ip.mat@free] <- NA

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

m2 <- mxModel(model="grm1", ip.mat, m.mat, cov.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=spec, ItemParam="itemParam",
                mean="mean", cov="cov", qpoints=31),
              mxFitFunctionML(),
	      mxComputeOnce('expectation', context='EM'))
middle <- mxRun(m2)
omxCheckCloseEnough(sum(middle@expectation@patternLikelihood), -9742.31, .1)
omxCheckCloseEnough(fivenum(middle@expectation@patternLikelihood),
                    c(-34.98313, -22.50933, -19.59691, -16.79153, -6.51683), .001)
omxCheckCloseEnough(sum(middle@expectation@em.expected), 7000, .01)
omxCheckCloseEnough(fivenum(middle@expectation@em.expected),
                    c(0, 0, 0.00451, 1.61901, 80.99209), .01)

testDeriv <- mxModel(m2,
	      mxComputeSequence(steps=list(
				 mxComputeOnce('expectation', context='EM'),
				 mxComputeOnce('fitfunction', fit=TRUE, gradient=TRUE, hessian=TRUE, ihessian=TRUE)
				 )))
testDeriv <- mxRun(testDeriv)
omxCheckCloseEnough(testDeriv@fitfunction@result, 9399.954, .01)
omxCheckCloseEnough(fivenum(testDeriv@output$gradient),
                    c(-14424.48407, -62.52714, -2.51876, 71.87544, 14651.19963), .01)
omxCheckCloseEnough(fivenum(testDeriv@output$hessian[testDeriv@output$hessian != 0]),
                    c(-1038404.94356, -20.82415, 0.06734, 53.01333, 1038503.00369 ), .01)
omxCheckCloseEnough(solve(testDeriv@output$hessian), testDeriv@output$ihessian, 1e-2)

plan <- mxComputeEM('expectation',
		    mxComputeNewtonRaphson(free.set='itemParam'),
		    mxComputeOnce('fitfunction', fit=TRUE, free.set=c("mean","cov")))

m2 <- mxModel(m2,
              mxExpectationBA81(
                ItemSpec=spec, ItemParam="itemParam",
                mean="mean", cov="cov",
                qpoints=31,
                scores="full"),
	      plan)
				  
m2 <- mxRun(m2)

#print(m2@matrices$itemParam@values)
#print(correct.mat)
omxCheckCloseEnough(m2@fitfunction@result, 13969.747, .01)
got <- cor(c(m2@matrices$itemParam@values[ip.mat@free]),
           c(correct.mat[ip.mat@free]))
omxCheckCloseEnough(got, .993, .01)
scores <- m2@expectation@scores.out
omxCheckCloseEnough(scores[1:5,1], c(0.81609, 0.74994, -0.83515, 0.79766, 0.16879), .001)
omxCheckCloseEnough(scores[1:5,2], c(0.43522, 0.44211, 0.4686, 0.43515, 0.3842), .001)
omxCheckCloseEnough(sum(abs(scores[,1] - ability) < 1*scores[,2])/500, .714, .01)
omxCheckCloseEnough(sum(abs(scores[,1] - ability) < 2*scores[,2])/500, .95, .01)
