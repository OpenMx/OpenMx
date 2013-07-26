# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
library(mvtnorm)

set.seed(15)
correct.deriv <- c(484.533329156793, 892.742069416003, 156.534557121965, 283.449559864665,  9.29526394590284)
correct.bifactor <- c(-48.70145, 280.13851, 23.39651, 109.05235, 31.38746, 209.36438,
                      172.0925, 124.71525, 164.80765, 131.46529)

numItems <- 20
numPeople <- 1000

items <- list()
correct <- list()
for (ix in 1:numItems) {
	items[[ix]] <- rpf.grm(factors=2)
	correct[[ix]] <- rpf.rparam(items[[ix]])
}
correct.mat <- simplify2array(correct)

maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))

true.mean <- c(.2,-.1)
true.cov <- matrix(c(.5, -.2, -.2, .5), nrow=2)
ability <- rmvnorm(numPeople, mean=true.mean, sigma=true.cov)
data <- rpf.sample(t(ability), items, correct.mat)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=correct.mat, free=FALSE)

eip.mat <- mxAlgebra(ItemParam, name="EItemParam")

m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=c(.5, .5), free=TRUE)
cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=matrix(c(1,.2,.2,1), nrow=2),
                    free=TRUE, labels=c("v1","c12","c12","v2"))

m1 <- mxModel(model="latent",
              ip.mat, m.mat, cov.mat, eip.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec=items,
                                EItemParam="EItemParam", ItemParam="ItemParam"),
              mxFitFunctionML(),
              mxComputeSequence(steps=list(
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction', gradient=TRUE))))

if (1) {
#  m1 <- mxOption(m1, "Number of Threads", 1)
#  m1 <- mxOption(m1, "No Sort Data", 'latent')
  m1 <- mxRun(m1, silent=TRUE)
  omxCheckCloseEnough(m1@expectation@empirical.mean, c(.168, .005), .01)
  omxCheckCloseEnough(m1@expectation@empirical.cov, matrix(c(.676, -.144, -.144, .682), nrow=2), .01)
  omxCheckCloseEnough(m1@output$gradient, correct.deriv, 1e-3)
}

if (1) {
  m1 <- mxModel(m1,
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation'),
                  mxComputeGradientDescent(useGradient=TRUE))))
#  m1 <- mxOption(m1, "Number of Threads", 1)
  m1 <- mxOption(m1, "Verify level", '-1')
  m1 <- mxOption(m1, "Function precision", '1.0E-5')
  m1 <- mxRun(m1, silent=TRUE)
  omxCheckCloseEnough(c(m1@matrices$mean@values), true.mean, .02)
  omxCheckCloseEnough(m1@matrices$cov@values, true.cov, .1)
}

objective1 <- function(x) {
	m.mat@values[1:2] <- x[1:2]
	cov.mat@values[1,1] <- x[3]
	cov.mat@values[1,2] <- x[4]
	cov.mat@values[2,1] <- x[4]
	cov.mat@values[2,2] <- x[5]
	m1 <- mxModel(m1, m.mat, cov.mat,
		      mxComputeSequence(steps=list(
					  mxComputeOnce("EItemParam"),
					  mxComputeOnce('expectation'),
					  mxComputeOnce('fitfunction'))))
	m1 <- mxRun(m1, silent=TRUE)
	got <- m1@output$minimum
#  print(got)
  got
}

if (0) {
  require(numDeriv)
  deriv <- grad(objective1, c(m.mat@values, cov.mat@values[1:2,1], cov.mat@values[2,2]))
  
  omxCheckCloseEnough(c(deriv), correct.deriv, 1e-3)
}

if (1) {
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=5, values=rep(.1,5), free=TRUE)
  cov.mat <- mxMatrix(name="cov", nrow=5, ncol=5, values=diag(seq(.81,1.41,length.out=5)),
                      free=diag(5)==1)
  
  design <- matrix(c(rep(1L,numItems),
                     as.integer(kronecker(2:5,rep(1,5)))), byrow=TRUE, ncol=numItems)
  m1 <- mxModel(m1, m.mat, cov.mat,
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=items,
                                  design=design,
                                  EItemParam="EItemParam", ItemParam="ItemParam"),
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation'),
                  mxComputeOnce('fitfunction', gradient=TRUE))))
#  m1 <- mxOption(m1, "Number of Threads", 1)
  m1 <- mxRun(m1, silent=TRUE)
  omxCheckCloseEnough(m1@output$gradient[1:5], correct.bifactor[1:5], 1e-3)  #means
  omxCheckCloseEnough(m1@output$gradient[6:10], correct.bifactor[6:10], 1e-2) #vars
}

objective2 <- function(x) {
  m.mat@values[1:5] <- x[1:5]
  cov.mat@values <- diag(x[6:10])
  m1 <- mxModel(m1, m.mat, cov.mat,
                mxComputeSequence(steps=list(
                  mxComputeOnce("EItemParam"),
                  mxComputeOnce('expectation'),
                  mxComputeOnce('fitfunction'))))
  m1 <- mxRun(m1, silent=TRUE)
  got <- m1@output$minimum
  #  print(got)
  got
}

if (0) {
  require(numDeriv)
  deriv <- grad(objective2, c(m.mat@values, diag(cov.mat@values)))
  cat(deparse(round(deriv, 5)))
  omxCheckCloseEnough(c(deriv), correct.deriv, 1e-3)
}
