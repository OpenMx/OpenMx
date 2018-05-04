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
