#
#   Copyright 2007-2026 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
# GREML (semi-)analytic gradients inside parallel-CI kid contexts.
#
# Numeric dV/dyhat derivatives use one child FitContext per OpenMP thread, but
# CI kids have no children of their own, and the derivative team historically
# numbered its threads with the nesting-absolute id (wrong work shares and
# childList indexes inside the CI team). Both must be innermost-team-relative
# and capped at the available child contexts (OpenMx #380: segfault at
# omxGREMLfitfunction.cpp:1970 under multiple threads).
#

library(OpenMx)
#mxOption(key='Number of Threads', value=4)
#library(testthat)

#context("GREML CIs with threads")

#test_that("GREML confidence interval matches serial with 4 threads", {
ompThreads = suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", unset = "4")))
if (is.na(ompThreads) || ompThreads < 1) ompThreads = 4L
nthreads = min(ompThreads, 4L)
if (nthreads < 2){ stop("SKIP") }

set.seed(476)
nobs = 1000L
grm = matrix(0, nobs, nobs)
grm[lower.tri(grm)] = runif(nobs * (nobs - 1) / 2, -0.025, 0.025)
grm = grm + t(grm)
diag(grm) = runif(nobs, 0.95, 1.05)
pheno = drop(chol(grm * 0.5) %*% rnorm(nobs)) + rnorm(nobs, sd = sqrt(0.5))
covar = rnorm(nobs)
dat = cbind(y = pheno, x = covar)

ciByThreads = list()
for (nthread in c(1L, nthreads)) {
	testmod = mxModel(
		"GREML_1GRM_1trait",
		mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE,
						 values = var(pheno) / 2, labels = "ve", lbound = 0.0001, name = "Ve"),
		mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE,
						 values = var(pheno) / 2, labels = "va", name = "Va"),
		mxMatrix("Iden", nrow = nobs, name = "I"),
		mxMatrix("Symm", nrow = nobs, free = FALSE, values = grm, name = "A"),
		mxAlgebra((A %x% Va) + (I %x% Ve), name = "V"),
		mxAlgebra(Va / (Va + Ve), name = "h2"),
		mxCI("h2"),
		mxData(observed = dat, type = "raw", sort = FALSE),
		mxExpectationGREML(V = "V", yvars = "y", Xvars = "x", addOnes = TRUE),
		mxFitFunctionGREML())
	#testmod$fitfunction$.parallelDerivScheme = 1L
	testmod = mxOption(testmod, "Number of Threads", nthread)
	testrun = mxRun(testmod, intervals = TRUE)
	omxCheckEquals(as.integer(testrun$output$status$code), 0L)
	gotCI = testrun$output$confidenceIntervals[1, c("lbound", "estimate", "ubound")]
	omxCheckTrue(all(is.finite(gotCI)))
	ciByThreads[[as.character(nthread)]] = gotCI
}
omxCheckEquals(ciByThreads[[as.character(nthreads)]], ciByThreads[["1"]], tolerance = 1e-4)
#})
