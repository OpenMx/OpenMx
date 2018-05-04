library(OpenMx)

# Rastrgin function (The global minimum is 0 with all parameters at 0.)
Rastrigin <- function(x) {
  sum(x^2 - 10 * cos(2 * pi  * x)) + 10 * length(x)
}

set.seed(1234)
dimension <- 8

m1 <- mxModel(
	"rastrigin",
	mxMatrix(type='Full', values=runif(dimension, min=-1e6, max=1e6),
		ncol=1, nrow=dimension, free=TRUE, name='x'),
	mxAlgebra(sum(x*x - 10 * cos(2 * pi * x)) + 10 * dimension, name="fit"),
	mxFitFunctionAlgebra('fit'),
	mxComputeGenSA(verbose=0L))

m1 <- mxRun(m1)

omxCheckCloseEnough(m1$output$fit, 0, 1e-2)
omxCheckCloseEnough(coef(m1), rep(0,dimension), 1e-2)

# ----

path <- paste0(tempdir(),'/','log.csv')

m2 <- mxModel(
	"rastrigin",
	mxMatrix(type='Full', values=runif(dimension, min=-1e6, max=1e6),
		ncol=1, nrow=dimension, free=TRUE, name='x'),
	mxAlgebra(sum(x*x - 10 * cos(2 * pi * x)) + 10 * dimension, name="fit"),
	mxFitFunctionAlgebra('fit'),
	mxComputeCheckpoint(path=path, append=FALSE))

m2 <- mxRun(m2)
m2 <- mxRun(m2)

log1 <- read.table(path, header=TRUE, sep='\t', check.names=FALSE)
omxCheckEquals(nrow(log1), 1)

m2$compute$append <- TRUE
m2$compute$header <- FALSE

m2 <- mxRun(m2)

log2 <- read.table(path, header=TRUE, sep='\t', check.names=FALSE)
omxCheckEquals(nrow(log2), 2)
