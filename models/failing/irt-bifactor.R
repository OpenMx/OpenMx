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
}
correct.mat <- simplify2array(correct)

maxParam <- max(vapply(items, rpf.numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i@outcomes, 0))

design <- matrix(c(rep(1,numItems),
		   rep(2,numItems/2), rep(3, numItems/2)), byrow=TRUE, nrow=2)

theta <- t(rmvnorm(numPersons, mean=rnorm(3, sd=.25)))
data <- rpf.sample(theta, items, correct, design)

spec <- mxMatrix(name="ItemSpec", nrow=6, ncol=numItems,
         values=sapply(items, function(m) slot(m,'spec')),
         free=FALSE, byrow=TRUE)

design <- mxMatrix(name="Design", nrow=maxDim, ncol=numItems,
		   values=design)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(1.414, 1, 0, 0, 1),
		   lbound=c(1e-6, 1e-6, -1e6, 0, 0),
		   free=c(rep(TRUE, 3), FALSE, FALSE))

#ip.mat@values[2,1] <- correct.mat[2,1]
#ip.mat@free[2,1] <- FALSE

m1 <- mxModel(model="bifactor",
          spec, design,
          ip.mat,
          mxData(observed=data, type="raw"),
          mxExpectationBA81(
	     ItemSpec="ItemSpec",
	     Design="Design",
	     ItemParam="ItemParam",
	    qpoints=29,
	    scores="full"),
          mxFitFunctionBA81()
)

m1 <- mxOption(m1, "Analytic Gradients", 'no')
if (1) {
	m1 <- mxOption(m1, "Analytic Gradients", 'yes')
	m1 <- mxOption(m1, "Verify level", '-1')
}
m1 <- mxOption(m1, "Function precision", '1.0E-5')
m1 <- mxOption(m1, "Calculate Hessian", "No")
m1 <- mxOption(m1, "Standard Errors", "No")

m1 <- mxRun(m1, silent=TRUE)
#print(correct.mat)
#print(m1@matrices$ItemParam@values)
got <- cor(c(m1@matrices$ItemParam@values), c(correct.mat))
omxCheckCloseEnough(got, .966, .01)
omxCheckCloseEnough(cor(c(m1@output$ability[1,]), c(theta[1,])), .774, .01)
omxCheckCloseEnough(cor(c(m1@output$ability[3,]), c(theta[2,])), .657, .01)
omxCheckCloseEnough(cor(c(m1@output$ability[5,]), c(theta[3,])), .683, .01)

omxCheckCloseEnough(sum(abs(m1@output$ability[3,] - theta[2,]) < 2*m1@output$ability[4,]), 919)
omxCheckCloseEnough(sum(abs(m1@output$ability[3,] - theta[2,]) < 3*m1@output$ability[4,]), 987)
