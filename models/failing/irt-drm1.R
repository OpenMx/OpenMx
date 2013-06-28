library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 30
i1 <- rpf.drm()
items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
  correct[[ix]][3] <- 0
}
correct.mat <- simplify2array(correct)

data <- rpf.sample(500, items, correct)

spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
         values=c(rep(mxLookupIRTItemModelID("drm1"), numItems),
		 rep(2, numItems),  # outcomes
		 rep(1, numItems)), # dim
         free=FALSE, byrow=TRUE)

ip.mat <- mxMatrix(name="itemParam", byrow=TRUE, nrow=3, ncol=numItems,
                   values=c(rep(1,numItems),
		     rep(0,numItems*2)),
                   free=c(
		     rep(TRUE,numItems),
		     rep(TRUE,numItems),
		     rep(FALSE,numItems)))

m2 <- mxModel(model="drm1", ip.mat, spec,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec="ItemSpec",
                ItemParam="itemParam",
		GHpoints=30),
              mxFitFunctionBA81())

m2 <- mxOption(m2, "Analytic Gradients", 'no')
if (1) {
	m2 <- mxOption(m2, "Analytic Gradients", 'yes')
	m2 <- mxOption(m2, "Verify level", '-1')
}
m2 <- mxOption(m2, "Calculate Hessian", "No")
m2 <- mxOption(m2, "Standard Errors", "No")
m2 <- mxRun(m2)

print(m2@matrices$itemParam@values)
print(correct.mat)
got <- cor(c(m2@matrices$itemParam@values),
           c(correct.mat))
print(got)
#omxCheckCloseEnough(got, .946, .01)

q()
