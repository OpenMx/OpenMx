#
# EM M-step-start cache invalidation for shared-state child contexts.
#
# The EM M-step runs under a child FitContext sharing the top state. It must
# still invalidate algebra/FIML caches at M-step start (without spawning
# grandchild teams). When that invalidation was skipped, this 2-class mixture
# fell into a +/-10% E/M 2-cycle from cycle 2, exhausted all 500 cycles
# (BLUE, status code 4) and landed at 535.38 instead of the optimum 539.1599.
#

library(OpenMx)
library(testthat)

context("EM M-step cache invalidation")

mxOption(key = "Number of Threads", value = 1)

test_that("mixture EM converges to the known optimum", {
	set.seed(190127)
	mu1 = -1.5
	mu2 = 1.5
	N = 200
	x = matrix(c(rnorm(N/2, mu1, 1), rnorm(N/2, mu2, 1)), ncol = 1, dimnames = list(NULL, "x"))
	data4mx = mxData(observed = x, type = "raw")
	class1 = mxModel("Class1",
		mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = -.5, name = "Mu"),
		mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = 4, name = "Sigma"),
		mxExpectationNormal(covariance = "Sigma", means = "Mu", dimnames = "x"),
		mxFitFunctionML(vector = TRUE))
	class2 = mxRename(class1, "Class2")
	class2$Mu$values[1,1] = .5
	mm = mxModel(
		"Mixture", data4mx, class1, class2,
		mxAlgebra((1 - Posteriors) * Class1.fitfunction, name = "PL1"),
		mxAlgebra(Posteriors * Class2.fitfunction, name = "PL2"),
		mxAlgebra(PL1 + PL2, name = "PL"),
		mxAlgebra(PL2 / PL, recompute = 'onDemand',
			initial = matrix(runif(N, .4, .6), nrow = N, ncol = 1), name = "Posteriors"),
		mxAlgebra(-2 * sum(log(PL)), name = "FF"),
		mxFitFunctionAlgebra(algebra = "FF"),
		mxComputeEM(
			estep = mxComputeOnce("Mixture.Posteriors"),
			mstep = mxComputeGradientDescent(fitfunction = "Mixture.fitfunction")))
	mmFit = mxRun(mm)
	expect_equal(as.integer(mmFit$output$status$code), 0L)
	expect_equal(mmFit$output$fit, 539.1599, tolerance = 0.01)
})
