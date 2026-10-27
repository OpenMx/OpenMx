#
# Selection-RAM finite-difference gradients under multithreading.
#
# RAM models with arrows=0 (Pearson selection) paths keep their free selection
# values in a selectionVector matrix. Gradient kids must read their own state's
# duplicate; when they aliased the parent's copy, perturbed evaluations all
# reused the unperturbed selection values, so gradients w.r.t. selection paths
# were exactly zero and the optimizer sat at its start values (OpenMx #380).
#

library(OpenMx)
library(testthat)

context("Pearson selection with threads")

test_that("selection RAM recovers truth with 4 threads", {
	ompThreads = suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "4")))
	if (is.na(ompThreads) || ompThreads < 1) ompThreads = 4L
	nthreads = min(ompThreads, 4L)
	if (nthreads < 2) skip("needs 2+ OpenMP threads")

	t3 <- mxModel(
		"t3",
		mxMatrix("Symm", 4, 4, FALSE, values = diag(4) + .2, name = "c1"),
		mxMatrix("Symm", 4, 4, free = c(F, T, rep(F, 8)),
			values = diag(4) + .2, labels = c(NA, 'p1', rep(NA, 8)), name = "c2"),
		mxAlgebra(mxPearsonSelCov(c1, c2), name = "c3"),
		mxMatrix("Full", 4, 1, values = .1, name = "m1"),
		mxAlgebra(mxPearsonSelMean(c1, c2, m1), name = "m2"),
		mxExpectationNormal(covariance = "c3", means = "m2"),
		mxFitFunctionML())

	trueSelCor = c(-.2, -.3)
	t3 = omxSetParameters(t3, "p1", values = trueSelCor[1])

	t3 = mxModel(t3,
		mxMatrix("Symm", 4, 4, free = c(F, F, T, rep(F, 7)),
			values = mxEval(c3, t3, compute = TRUE),
			labels = c(NA, NA, 'p2', rep(NA, 7)), name = "c4"),
		mxAlgebra(mxPearsonSelCov(c3, c4), name = "c5"),
		mxAlgebra(mxPearsonSelMean(c3, c4, m2), name = "m3"),
		mxExpectationNormal(covariance = "c5", means = "m3"))

	t3 = omxSetParameters(t3, "p2", values = trueSelCor[2])

	set.seed(1)
	dat = mxGenerateData(t3, nrows = 800)

	t4 = mxModel(
		"t4", type = "RAM",
		manifestVars = paste0('V', 1:4),
		mxPath(paste0('V', 1:4), arrows = 2, values = 1.2, free = FALSE),
		mxPath(paste0('V', 1:4), arrows = 2,
			connect = "unique.bivariate", values = .2, free = FALSE),
		mxPath('one', paste0('V', 1:4), values = .1, free = FALSE),
		mxPath('V1', 'V2', arrows = 0, values = .2),
		mxPath('V1', 'V3', arrows = 0, values = .2, step = 2),
		mxData(dat, 'raw'))
	t4 = mxOption(t4, "Number of Threads", nthreads)

	t4 = mxRun(t4)

	expect_equivalent(coef(t4), trueSelCor, .03)
})
