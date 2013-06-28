# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
#library(mvtnorm)

set.seed(7)

numItems <- 5
numPeople <- 500
maxDim <- 2

items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
	items[[ix]] <- rpf.drm(factors=min(ix,maxDim),
                               multidimensional=TRUE)
	correct[[ix]] <- rpf.rparam(items[[ix]])
	if (ix>1) correct[[ix]][[4]] <- 0   # no guessing, for now
}
correct[[1]][[5]] <- 1   # make all vectors the same length
correct.mat <- simplify2array(correct)
correct.mat[5,] <- 1
correct.mat[4,1] <- 1
correct.mat[3,1] <- 0

maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))
maxOutcomes <- max(vapply(items, function(i) i@outcomes, 0))

design <- matrix(c(1, 1,1,1,2,
                  NA,2,2,2,1), byrow=TRUE, nrow=2)

ability <- array(rnorm(numPeople * 2), dim=c(2, numPeople))
#cov <- matrix(c(1, .68, .68, 1), nrow=2)
#ability <- rmvnorm(numPeople, sigma=cov)
data <- rpf.sample(ability, items, correct.mat, design)

spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
         values=sapply(items, function(m) slot(m,'spec')),
         free=FALSE, byrow=TRUE)

design <- mxMatrix(name="Design", nrow=maxDim, ncol=numItems,
                  values=design)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=c(1, 1.4, 0, 0, 1),
                  lbound=c(1e-6, -1e6, 0,  0, 0,
                    rep(c(1e-6, 1e-6, -1e6, 0, 0), 4)),
         free=c(TRUE, TRUE, FALSE, FALSE, FALSE,
          rep(c(TRUE, TRUE, TRUE, FALSE, FALSE), 4)))
ip.mat@values[4,1] <- 1

m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=diag(2), free=FALSE)

m1 <- mxModel(model="2dim",
          spec, design,
          ip.mat, m.mat, cov.mat,
          mxData(observed=data, type="raw"),
          mxExpectationBA81(mean="mean", cov="cov",
	     ItemSpec="ItemSpec",
	     Design="Design",
	     ItemParam="ItemParam",
	    qpoints=29,
	    scores="full"),
          mxFitFunctionBA81())

m1 <- mxOption(m1, "Analytic Gradients", 'no')
if (1) {
	m1 <- mxOption(m1, "Analytic Gradients", 'yes')
	m1 <- mxOption(m1, "Verify level", '-1')
}
m1 <- mxOption(m1, "Function precision", '1.0E-5')
m1 <- mxOption(m1, "Calculate Hessian", "No")
m1 <- mxOption(m1, "Standard Errors", "No")

m1 <- mxRun(m1, silent=TRUE)

					#print(m1@matrices$ItemParam@values)
					#print(correct.mat)
# sometimes found as low as .89, maybe solution is unstable
omxCheckCloseEnough(cor(c(m1@matrices$ItemParam@values),
			c(correct.mat)), .91, .04)
max.se <- max(m1@output$ability[c(2,4),])
omxCheckCloseEnough(sum(abs(m1@output$ability[c(1,3),] - ability) < max.se) / (numPeople*2), .797, .01)
omxCheckCloseEnough(.672, cor(c(m1@output$ability[c(1,3),]), c(ability)), .01)
