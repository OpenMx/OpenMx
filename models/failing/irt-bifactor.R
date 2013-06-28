# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)

set.seed(5)

numItems <- 5
numPersons <- 1000
maxDim <- 2

items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
	items[[ix]] <- rpf.drm(dimensions=maxDim)
	correct[[ix]] <- rpf.rparam(items[[ix]])
	correct[[ix]][[4]] <- 0   # no guessing, for now
}
correct.mat <- simplify2array(correct)

maxParam <- max(vapply(items, function(i) i@numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i@numOutcomes, 0))

design <- matrix(c(rep(1,numItems),
		   2,2,2,3,3), byrow=TRUE, nrow=2)

data <- rpf.sample(numPersons, items, correct, design)

spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
         values=c(rep(mxLookupIRTItemModelID("drm"), numItems),
		 rep(2, numItems),
		 rep(2, numItems)),
         free=FALSE, byrow=TRUE)

design <- mxMatrix(name="Design", nrow=maxDim, ncol=numItems,
		   values=design)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(
		     rep(1.414,numItems),
		     rep(1,numItems),
		     rep(0,numItems),
		     rep(0,numItems)),
         free=c(rep(TRUE, numItems*3),
		 rep(FALSE, numItems*1)),
         byrow=TRUE)

#ip.mat@values[2,1] <- correct.mat[2,1]
#ip.mat@free[2,1] <- FALSE

m1 <- mxModel(model="bifactor",
          spec, design,
          ip.mat,
          mxData(observed=data, type="raw"),
          mxExpectationBA81(
	     ItemSpec="ItemSpec",
	     Design="Design",
	     ItemParam="ItemParam"),
          mxFitFunctionBA81()
)

m1 <- mxOption(m1, "Analytic Gradients", 'no')
if (0) {
	m1 <- mxOption(m1, "Analytic Gradients", 'yes')
	m1 <- mxOption(m1, "Verify level", '-1')
}
m1 <- mxOption(m1, "Calculate Hessian", "No")
m1 <- mxOption(m1, "Standard Errors", "No")

m1 <- mxRun(m1, silent=TRUE)
print(correct.mat)
print(m1@matrices$ItemParam@values)
got <- cor(c(m1@matrices$ItemParam@values), c(correct.mat))
omxCheckCloseEnough(got, .957, .01) # caution, not converged yet
