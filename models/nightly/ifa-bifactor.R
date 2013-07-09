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

design <- matrix(c(rep(1,numItems),
		   rep(2,numItems/2), rep(3, numItems/2)), byrow=TRUE, nrow=2)

theta <- t(rmvnorm(numPersons, mean=rnorm(3, sd=.25)))
data <- rpf.sample(theta, items, correct, design)

spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
         values=sapply(items, function(m) slot(m,'spec')),
         free=FALSE, byrow=TRUE)

design <- mxMatrix(name="Design", nrow=maxDim, ncol=numItems,
		   values=design)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(1.414, 1, 0, 0, 1),
		   lbound=c(1e-6, 1e-6, -1e6, 0, 0),
		   free=c(rep(TRUE, 3), FALSE, FALSE))
ip.mat@free.group <- 'param'

#ip.mat@values[2,1] <- correct.mat[2,1]
#ip.mat@free[2,1] <- FALSE

eip.mat <- mxAlgebra(ItemParam, name="EItemParam")

m.mat <- mxMatrix(name="mean", nrow=1, ncol=3, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=3, ncol=3, values=diag(3), free=FALSE)

m1 <- mxModel(model="bifactor",
          spec, design,
          ip.mat, m.mat, cov.mat, eip.mat,
          mxData(observed=data, type="raw"),
          mxExpectationBA81(mean="mean", cov="cov",
	     ItemSpec="ItemSpec",
	     Design="Design",
	     EItemParam="EItemParam",
	    qpoints=29,
	    scores="full"),
          mxFitFunctionBA81(ItemParam="ItemParam"),
              mxComputeIterate(steps=list(
                mxComputeOnce("EItemParam"),
                mxComputeOnce('expectation', context='E'),
			   mxComputeNewtonRaphson(free.group='param'),
#                mxComputeGradientDescent(free.group='param'),
                mxComputeOnce('expectation', context='M'),
                mxComputeOnce('fitfunction'))))

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
