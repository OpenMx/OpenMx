# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)

set.seed(7)

numItems <- 5
numPersons <- 500
maxDim <- 2

items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
	items[[ix]] <- rpf.drm(dimensions=min(ix,maxDim),
                               multidimensional=TRUE)
	correct[[ix]] <- rpf.rparam(items[[ix]])
	correct[[ix]][[4]] <- 0   # no guessing, for now
}
correct[[1]][[3]] <- 0   # no guessing, for now
correct.mat <- simplify2array(correct)

maxParam <- max(vapply(items, function(i) i@numParam, 0))
maxOutcomes <- max(vapply(items, function(i) i@numOutcomes, 0))

design <- matrix(c(1, 1,1,1,2,
		   NA,2,2,2,1), byrow=TRUE, nrow=2)

data <- rpf.sample(numPersons, items, correct, design)

colnames(data) <- paste("item", 1:dim(data)[2], sep='')

spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
         values=c(rep(mxLookupIRTItemModelID("drm"), numItems),
		 rep(2, numItems),
		 c(1,2,2,2,2)),
         free=FALSE, byrow=TRUE)

design <- mxMatrix(name="Design", nrow=maxDim, ncol=numItems,
		   values=design)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(
		     rep(1,numItems),
		     rep(1,numItems),
		     rep(0,numItems),
		     rep(0,numItems)),
         free=c(rep(TRUE, numItems),
	        rep(TRUE, numItems),
	        FALSE, rep(TRUE, numItems-1),
		 rep(FALSE, numItems)),
         byrow=TRUE)

m1 <- mxModel(model="2dim",
          spec, design,
          ip.mat,
          mxMatrix(name="A", values=c(1,1,0,0), nrow=1, ncol=maxParam, free=FALSE),
          mxAlgebra(name="prior",
                    sum(dlnorm(omxSelectRows(ItemParam, A), 0, .5, 1))),
          mxData(observed=data, type="raw"),
          mxExpectationBA81(
	     ItemSpec="ItemSpec",
	     Design="Design",
	     ItemParam="ItemParam",
	     ItemPrior="prior"),
          mxFitFunctionBA81())

m1 <- mxRun(m1, silent=TRUE)

omxCheckCloseEnough(cor(c(m1@matrices$ItemParam@values),
			c(correct.mat)), .93, .01)

################
if (0){
myrpf <- function(iparam, spec, where) {
	ret <- array(dim=c(maxOutcomes, numItems))
	for (ix in 1:numItems) {
		i <- items[[ix]]
		param <- iparam[1:i@numParam,ix]
		out <- rpf.logprob(i, param, t(where[1:i@dimensions,ix]))
		ret[1:i@numOutcomes,ix] <- out
	}
	ret
}

m2 <- mxModel(name="R", m1,
	      mxExpectationBA81(
		ItemSpec="ItemSpec",
		Design="Design",
		ItemParam="ItemParam",
		ItemPrior="prior",
		RPF=myrpf))

m2 <- mxRun(m2, silent=TRUE) # 23min, ugh
omxCheckCloseEnough(cor(c(m2@matrices$ItemParam@values),
			c(correct.mat)), .93, .01)
}
