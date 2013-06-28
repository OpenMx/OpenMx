# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)

set.seed(7)

numItems <- 5
numPeople <- 500
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

ability <- array(rnorm(numPeople * 2), dim=c(numPeople, 2))
data <- rpf.sample(ability, items, correct, design)

spec <- mxMatrix(name="ItemSpec", nrow=6, ncol=numItems,
         values=sapply(items, function(m) slot(m,'spec')),
         free=FALSE, byrow=TRUE)

design <- mxMatrix(name="Design", nrow=maxDim, ncol=numItems,
		   values=design)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(
		     rep(1,numItems),
		     rep(1.4,numItems),
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
          mxData(observed=data, type="raw"),
          mxExpectationBA81(
	     ItemSpec="ItemSpec",
	     Design="Design",
	     ItemParam="ItemParam",
	    GHpoints=21),
          mxFitFunctionBA81())

m1 <- mxOption(m1, "Analytic Gradients", 'no')
if (1) {
	m1 <- mxOption(m1, "Analytic Gradients", 'yes')
	m1 <- mxOption(m1, "Verify level", '-1')
}
m1 <- mxOption(m1, "Function precision", '1.0E-5')
m1 <- mxOption(m1, "Calculate Hessian", "No")
m1 <- mxOption(m1, "Standard Errors", "No")

if (1) {

	m1 <- mxRun(m1, silent=TRUE)

	omxCheckCloseEnough(cor(c(m1@matrices$ItemParam@values),
				c(correct.mat)), .936, .01)
	max.se <- max(m1@output$ability[,c(2,4)])
	omxCheckCloseEnough(m1@output$ability[,c(1,3)], ability, max.se*3)
	omxCheckCloseEnough(.587, cor(c(m1@output$ability[,c(1,3)]), c(ability)), .01)

} else { ################

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
				c(correct.mat)), .934, .01)
}

q()
