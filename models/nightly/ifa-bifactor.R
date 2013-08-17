# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
library(mvtnorm)

set.seed(5)

numItems <- 20
numPersons <- 1000
maxDim <- 2

items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
	items[[ix]] <- rpf.drm(factors=maxDim)
	correct[[ix]] <- rpf.rparam(items[[ix]])
	correct[[ix]][[4]] <- 0   # no guessing, for now
	correct[[ix]][[5]] <- 1   # upper cound
}
correct.mat <- simplify2array(correct)

maxParam <- max(vapply(items, rpf.numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i@outcomes, 0))

design <- matrix(c(rep(1L,numItems),
		   rep(2L,numItems/2), rep(3L, numItems/2)), byrow=TRUE, nrow=2)

theta <- t(rmvnorm(numPersons, mean=rnorm(3, sd=.25)))
data <- rpf.sample(theta, items, correct, design)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(1.414, 1, 0, 0, 1),
		   free=c(rep(TRUE, 3), FALSE, FALSE))

#ip.mat@values[2,1] <- correct.mat[2,1]
#ip.mat@free[2,1] <- FALSE

m.mat <- mxMatrix(name="mean", nrow=1, ncol=3, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=3, ncol=3, values=diag(3), free=FALSE)

m1 <- mxModel(model="bifactor",
          ip.mat, m.mat, cov.mat,
          mxData(observed=data, type="raw"),
          mxExpectationBA81(mean="mean", cov="cov",
	     ItemSpec=items,
	     design=design,
	     ItemParam="ItemParam",
	    qpoints=29),
	      mxFitFunctionML(),
	      mxComputeOnce('expectation', context='EM'))
m1 <- mxRun(m1)

omxCheckCloseEnough(sum(m1@expectation@patternLikelihood), -12629.4, .1)
omxCheckCloseEnough(fivenum(m1@expectation@patternLikelihood),
                    c(-15.7575854, -14.9684791, -14.0992631, -12.3467773, -3.5902924 ), 1e-4)
omxCheckCloseEnough(sum(m1@expectation@em.expected), 20000, 1)
omxCheckCloseEnough(fivenum(m1@expectation@em.expected),
                    c(0, 0, 1.8e-06, 0.0034365, 43.2895967), 1e-4)

m1 <- mxModel(m1,
              mxComputeIterate(steps=list(
                mxComputeOnce('expectation', context='EM'),
                mxComputeOnce('fitfunction', fit=TRUE, gradient=TRUE, hessian=TRUE)
              )))
m1 <- mxRun(m1)
omxCheckCloseEnough(m1@fitfunction@result, 11850.68, .01)
omxCheckCloseEnough(fivenum(m1@output$gradient), c(-369.32879, -14.47296, 13.1165, 50.07066, 323.04627 ), .01)
omxCheckCloseEnough(fivenum(m1@output$hessian[m1@output$hessian != 0]),
                    c(-53.666201, -7.6857353, -6.0121325, 89.8735155, 192.6600613 ), 1e-4)

m1 <- mxModel(m1,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov", cache=FALSE,
                                ItemSpec=items,
                                design=design,
                                ItemParam="ItemParam",
                                qpoints=29, scores="full"),
              mxComputeIterate(steps=list(
                mxComputeOnce('expectation', context='EM'),
			   mxComputeNewtonRaphson(free.set='ItemParam'),
#                mxComputeGradientDescent(free.set='ItemParam'),
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction', free.set=c("mean","cov"), fit=TRUE)
				 )))

	m1 <- mxOption(m1, "Analytic Gradients", 'Yes')
	m1 <- mxOption(m1, "Verify level", '-1')
m1 <- mxOption(m1, "Function precision", '1.0E-5')

m1 <- mxRun(m1, silent=TRUE)
#print(correct.mat)
#print(m1@matrices$ItemParam@values)
got <- cor(c(m1@matrices$ItemParam@values), c(correct.mat))
omxCheckCloseEnough(got, .977, .01)
scores <- m1@expectation@scores.out
omxCheckCloseEnough(cor(c(scores[,1]), c(theta[1,])), .758, .01)
omxCheckCloseEnough(cor(c(scores[,2]), c(theta[2,])), .781, .01)
omxCheckCloseEnough(cor(c(scores[,3]), c(theta[3,])), .679, .01)

omxCheckCloseEnough(sum(abs(scores[,2] - theta[2,]) < 2*scores[,5]), 933, 5)
omxCheckCloseEnough(sum(abs(scores[,2] - theta[2,]) < 3*scores[,5]), 1000, 2)
