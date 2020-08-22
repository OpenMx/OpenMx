# Where to find news

OpenMx developers, being lazy and incorrigible, often forget to update the NEWS file. To learn about new and exciting features, please visit https://openmx.ssri.psu.edu/

# OpenMx 2.18
* October 2020 (R 4.0.2)
* FIXED: Bug in `mxGenerateData` where use.miss was ignored if nrowsProportion was set.


# OpenMx 2.17
* February 2020 (R ???)
 
# OpenMx 2.14.11
* October 2019 (R 3.6.1)
* IMPROVED: `mxRun` informs user about any unrecognised parameters (... parameters).
* FIXED: Bug affecting RAM models: `cov` data was not reordered correctly to match model.
* FIXED: Bug prevented estimation of factor score when there was more than one latent factor (v2.14.11-87-g6af1bf6fe).

# OpenMx 2.13.3
* September 2019 (R 3.6.1)
* FIXED: Bug `mxFactor`where RAM regression scores did not always include the means.
* FIXED: Bug `mxPowerSearch` was ignoring the user N. 
* FIXED: Bug `mxPower` where we were not respecting requested p-values.
* IMPROVED: `mxGenerateData` returns data of the type in the model (when `returnModel = TRUE``)
* IMPROVED: `mxPower` Clarify N reported by mxPower is total N, not average/group.
* IMPROVED: `mxPowerSearch` supports un-equal N in multi-group models, as specified by user's trueModel.

# OpenMx 2.13.2 
* May 2019 (R 3.6.0)
* NEW!: SEs for models with constraints! (used not to be calculated. Big win for twin models)
* NEW: `vcov` works with models that contain MxConstraints (so does `m1$vcov')
* EXPERIMENTAL: `mxAlgebra` gives users the option of populating an MxAlgebra with initial values, and of only recomputing the MxAlgebra when specifically called-for
* HELP: `mxRefModels` now raises a warning if it is used on an MxModel that contains definition variables, or is a multilevel model.
* IMPROVED: WLS summary statistics now taking advantage of multiple CPU cores: 40x faster :-)
* IMPROVED: `summary` now reports a p-value of 1 if the chi-square test statistic has zero degrees-of-freedom.
* NEW: `omxReadGRMBin` loads a genomic-relatedness matrix into R's workspace from a GCTA-format binary file.
* IMPROVED: `omxSetParameters` warns if you don't ask to do anything.
* NEW: `mxJiggle`: emulate the effect of keyword JIGGLE from classic Mx, or (default) function as a wrapper to the pre-existing imxJiggle(). Let the jiggling commence.
* IMPROVED: `mxData` allows non-positive-definite observed covariance matrix if type = "acov".
* IMPROVED: `mxRun` progress printing is now much easier to read
* IMPROVED: `mxTryHard` progress reporting much improved!
* IMPROVED: `mxOption` model  now conveniently defaults to NULL.
* IMPROVED: `mxOption` has a `reset` argument, to reset all options to their on-load defaults.
* EXPERIMENTAL: `mxComputeLoadData` supports high-throughput data analyses. A paper describing how to use it to analyze molecular genetic data is in preparation.
* EXPERIMENTAL: `mxData` now has a new argument, algebra. It is an experimental feature that is only useful in conjunction with mxComputeLoadData.
* IMPROVED: `mxCheckIdentification` works with models containing constraints
* IMPROVED: `mxFactorScores` is moved to the backend resulting in a substantial performance improvement!
* IMPROVED: `mxPath` error messages are even more helpful :-)
* IMPROVED: `mxModel` error when path added to model without type = "RAM"
* IMPROVED: Github and CRAN README.md
* CHANGED: CXX14 (not CXX11)


# OpenMx 2.12.2
* Feb 2019 (R 3.5.2)
* IMPROVED: parallel processing improved by addition of ConcurrentQueue.
* IMPROVED: WLS moved to the backend: 10x faster!
* NEW:  A new optimizer, which uses generalized simulated annealing, has been implemented.
* NEW:  `omxAkaikeWeights` and `mxModelAverage`, for information-theoretic model-averaging and multimodel inference.
* NEW:  `mxPearsonSelCov` and `mxPearsonSelMean` implement the Pearson-Aitken selection formulae. Both functions are usable in MxAlgebras.
* NEW:  `mxComputeLoadMatrix`: placed into a custom compute plan will load a CSV file directly into the backend.
* NEW:  `mxFitFunctionWLS` WLSM and WLSMV fit statistics.
* IMPROVED:  "CI details" table (in verbose summary() output) reordered, to make the table more readable.
* IMPROVED:  `mxStandardizeRAMpaths` now reports elements of the 'M' matrix, re-scaled to standard-deviation units.
* IMPROVED:  `mxAutoStart` can now be used with diagonally weighted least squares.
* IMPROVED:  `mxGenerateData` now compatible with models that depend on objects in other models.
* IMPROVED:  `mxOption`  "Max minutes" sets a maximum allowed backend time (default = 0, meaning no limit, i.e. Inf).
* IMPROVED:  `mxMatrix` now partially matches the value of its type argument. For instance, type="Ze" is now equivalent to type="Zero".
* IMPROVED:  `omxManifestModelByParameterJacobian` Jacobian  has dimnames, which make it easier to read.
* IMPROVED:  `SLSQP` uses multiple threads
* IMPROVED:  custom compute plan support checkpointing steps.
* IMPROVED:  Eigen linear-algebra library now multi-threaded.
* IMPROVED:  `mxExpectationHiddenMarkov` and `mxExpectationMixture` allow scaling of zero.
* IMPROVED:  `mxConstraints` that depend upon definition variables now throw a warning at runtime.
* IMPROVED:  `mxPower` more detailed output.
* IMPROVED:  `mxFitFunctionML`with rowDiagnostics=TRUE includes per-row squared Mahalanobis distance.

# OpenMx 2.12.1
* Jan 20 2019 (R 3.5.2)
* Bug fix to definition variable handling in multilevel models (v2.11.5-2-g7ef03e3fe)

# OpenMx 2.11.4
* September 24 2018 (R 3.5.1)
* PARTY TIME!:  Appears to be passing compiling for MacOS on CRAN !! 
* NEW: `mxModelAverage` function do compute parameter estimates that reflect the values found in a range of models that contain the parameter.
* IMPROVED: `mxTryHard` has compact and self-erasing progress report
* IMPROVED: `mxSE` is now **MUCH** faster - moved to the backend.
* IMPROVED: CIs on the RMSEA statistic for models that fit *very* badly.
* IMPROVED: Mahalanobis distance to ML models (Resolves #92)
* IMPROVED: `logLik.MxModel` can take a list of models.
* IMPROVED: `omxGetParameters` handles labels of the form model.mat[row, col]
* IMPROVED: Fit value is smart enough to not report logLik() for WLS (only reports AIC & BIC with -2lnL fit units)
* MODIFIED: `mxBrownie` now supports vegans.
* PREVIEW: `mxCompareMatrix` (note: this doesn't compare matrices, it compares models, and outputs a matrix of comparisons). Comments welcome!

Some other functions and changes that might interest you

* OTHER: `mxComputeManifestByParJacobian`

# OpenMx 2.9.9
* R 3.5.0) Summer 2018
* No change for users

# OpenMx 2.9.6
* R 3.4.4
* BUGFIX: Starting in version 2.7.0 (released January 10th, 2017), and ONLY in joint FIML analysis of both ordinal-threshold and continuous variables: if there was at least one row of the dataset in which all of the continuous variables had missing scores, then OpenMx could incorrectly evaluate the fit function, which would result in OpenMx silently returning numerically incorrect results. If you have run a joint ordinal/continuous FIML analysis with OpenMx versions 2.7.x thru 2.8.x, check your dataset to see if there are any rows in which all the continuous endogenous variables are missing (and at least one endogenous binary or ordinal variable is non-missing). If so, re-run your analysis using OpenMx v2.9. We apologize for any inconvenience this serious, but now-repaired, bug may have caused.
* IMPROVED:  CSOLNP. This bug would cause CSOLNP to enter an endless loop when optimizing an MxModel containing MxConstraints.
* IMPROVED:  Four memory leaks have been closed: one in OpenMx's interface to SLSQP, one in OpenMx's Nelder-Mead implementation, one in the MxComputeNumericDeriv step, and one in the FIML fitfunction backend for state-space models.
* BUGFIX: OpenMx's interface to SLSQP has been repaired. This bug caused SLSQP to always use the value of the 'Optimality tolerance' mxOption when checking the optimality conditions of a solution subject to MxConstraints. Now, as designed, if a user creating a custom compute plan passes a value for the optimality tolerance to mxComputeGradientDescent() that differs from the option, that different value controls.
* BUGFIX: All the known issues in the v2.8.3 release announcement have been resolved.
* NEW: Two new functions for power analysis, mxPower() and mxPowerSearch(), have been implemented.
* NEW: CSOLNP can now use analytic gradients and analytic constraint Jacobians, and can now numerically calculate gradients using the central-differences approximation.
* IMPROVED: `mxGenerateData` and parametric bootstrapping (including the bootstrap LRT) are now compatible with MxExpectationMixture.
* IMPROVED: Previously, nonparametric bootstrapping was incompatible with the FIML fitfunction when running on a single thread, or when using argument rowwiseParallel=FALSE to mxFitFunctionML().
* IMPROVED: When the 'maxOrdinalPerBlock' mxOption is too small, the FIML fitfunction backend no longer silently ignores nonzero covariance elements in order to coerce the dimension of a block of correlated ordinal variables to 'maxOrdinalPerBlock' or smaller.
* IMPROVED: When used on an MxModel that uses GREML expectation, mxGetExpected() now returns a vector of fitted values ("yhats" from feasible-generalized-least-squares regression) for "means", and now filters out rows and columns in the "covariance" matrix that correspond to missing observations.
* IMPROVED: Bias-corrected bootstrap-quantile confidence intervals will no longer sporadically have a lower limit of -Inf when the number of replications is too small for the desired coverage probability.
* IMPROVED: The print method for MxPaths has been improved.
* IMPROVED: `mxStandardizeRAMpaths` now throws a helpful error message if it detects that the MxModel contains free parameters the labels of which do not appear in the dimnames of the Hessian matrix.
* IMPROVED: `mxConstraint` documentation has been clarified regarding analytic constraint Jacobians for inequality constraints, and how MxConstraints interact with definition variables. Additionally, creating a constraint function that depends upon definition variables now raises a warning at runtime.
* IMPROVED: A potentially confusing typo in an MxExpectationLISREL warning message has been corrected.
* FIXED: It is now possible to build an NPSOL-enabled OpenMx binary with gcc version 7.x.
* NEW:  `omxModelDeleteData` which deletes data from an MxModel and all of its submodels recursively.
* IMPROVED: The RAM expectation no longer throws a runtime error if all of the diagonal elements of an 'S' matrix are zero (as that can be a reasonable scenario when conducting likelihood-ratio tests with multilevel models).
* IMPROVED: MxModel summary() output now includes the sample-size corrected AIC.
* IMPROVED: With argument details=TRUE, imxRobustSE() includes additional details in its output.

# OpenMx 2.8.3 (November 20, 2017)
* IMPROVED: `mxEval` can use scalar multiplication, division, and powering of matrices.
* IMPROVED: Ordinal data fit should often be better, and fewer status Reds.
* FIXED: mxBootstrap bug where replications reported incorrect optimizer status codes
* IMPROVED: Parametric bootstrapping compatible with multi-group models.
* IMPROVED: error checking in `mxCompare`, `mxCompareMatrix`, `confint`, `vcov`
	* Catch case where MxModel hasn't been run or has been modified since it was run.
* IMPROVED: `omxParallelCI` and `omxRunCI` accept new 'optimizer' argument for calculating CIs
* IMPROVED: All-ordinal models with complete data can return saturated multinomial model from `mxRefModels`.
* FIXED: Bug in `mxGetExpected`
* HELPFUL: OpenMx shows how to set the number of processor threads.
* NOTE: `mxBootstrap` may incorrectly report optimizer status codes for models using Nelder-Mead.
* NOTE: large number of replications with `mxBootstrap` may cause protect-stack errors.
* NEW: Add mxParametricBootstrap
* NEW: Add simulate S3 alias for mxGenerateData
* NEW: Add MxExpectationMixture interface.
* NEW: mxBootstrap, mxComputeBootstrap, mxBootstrapEval, mxBootstrapEvalByName
* NEW: mxEvaluateOnGrid algebra op (useful for quadrature integration)
* NEW: mxRobustLog algebra op (recommended for mixtures)
* IMPROVED: Nelder-Mead "GDsearch" method for handling equality constraints.
* IMPROVED: mxAutoStart handles models which refer to parent models.
* IMPROVED: WLS improvements
	* Handling models with constraints
	* sat models for mxFitFunctionMultigroup fit
* IMPROVED: mxComputeEM is more generic
* CHANGE: mxBootstrap no longer accepts a plan= parameter.
* FIXED: Oakes standard errors (ComputeEM) were broken since Feb 2016.
* FIXED: Bug in which NAs were allowed in integer definition variables.
* FIXED: Rampart sufficient statistic optimization no longer interferes with zero variance predictors.

# OpenMx 2.7.9 (March 22, 2017)
* NEW: mxAutoStart Get automatic starting values.
* NEW: "Auto" mxOption values.
* NEW: `omxNudgeZeroStarts` helper function
* NEW: Hidden Markov models with `mxExpectationHiddenMarkov`
* NEW: "Internal" warm starts for NPSOL
* NEW: Correct documentation of how NPSOL uses feasibility tolerance
* NEW: Nelder-Mead
* IMPROVED: mxEval() works with square brackets.
* IMPROVED: Summary table is easier to read and more informative
* IMPROVED: Rampart uses much less memory on large datasets.
* IMPROVED: 'one' is now in the list of reserved names (clashes with RAM's 'one' for mean)
* FIXED: mxKalmanScores are now correct for continuous time state space models.
* FIXED: Bugs in ref models with cov data.
* FIXED: Patched 3 bugs in Rampart's sufficient statistic optimization.
	* If using 2.7.x (older), set expectation$.useSufficientSets = FALSE to avoid this bug.

# OpenMx 2.7.x: Lots of NEWs
* NEW: mxSE Calculate standard errors for arbitrary named entities (free parameters, matrices, algebras) and expressions.
* NEW: mxRun Now gives feedback about optimization progress.
* NEW: Analytic constraint Jacobians can now be provided to MxConstraints (presently with NPSOL only).
* NEW: Constraint diagnostics exported to frontend (NPSOL and SLSQP).
* NEW: `omxDefaultComputePlan`
* CHANGE: CSOLNP is now the default optimizer.
	* Change this with mxOption(NULL, "Default optimizer", "CSOLNP|NPSOL|SLSQP")
	* NPSOL is available from the custom build at http://openmx.ssri.psu.edu/installing-openmx 
* SPEEDUP: In many cases, continuous data now evaluated as fast as covariance data instead of being much slower (SEM models including RAM and LISREL).
* SPEEDUP: F (means) matrix is optimised more efficiently.
* SPEEDUP: Evaluation of the GREML fitfunction's derivatives is now faster.
* IMPROVED: Error reporting is improved!
* IMPROVED: mxGenerateData can now take a data.frame instead of a model (returns data based on a saturated multivariate Gaussian model).
* IMPROVED: mxGenerateData now succeeds for joint ordinal and continuous data.
* INCOMPATIBLE CHANGE: mxGenerateData by default now approximates the missingness pattern of the original data. This can be turned off with use.miss=FALSE to get the previous behavior.
* INCOMPATIBLE CHANGE: mxFactorScores now asks user to specify the minimum number of values not NA in order to generate a score for a row.
* IMPROVED: Confidence interval diagnostics. See summary(..., verbose = TRUE).
* IMPROVED: Dynamically balance work between threads using empirical elapsed times when evaluating raw data in a row-wise parallel manner.
* IMPROVED: Continuous time state space models now allow non-invertible ("drift" or "dynamics") A matrices.
* IMPROVED: SLSQP now ignores inactive equality constraints and correctly report when inequality constraints cannot be satisfied.
* IMPROVED: mxCI() Wu & Neale (2012) adjustment for parameters with upper or lower bounds. 
	* Use mxCI(..., boundAdj = FALSE) to disable the adjustment.
* FIXED: mxRefModels Now handles models with one data variable.
* FIXED: mxTryHard Now no longer computes the Hessian and standard errors when there are MxConstraints in the model, which makes its behavior consistent with mxRun().
* NEW: Functions now usable in MxAlgebras: dchisq(), pchisq(), dbinom(), pbinom(), dcauchy(), pcauchy().
* IMPROVED: better man pages for mxOption() and mxComputeGradientDescent().
* FIXED: bug with omxSetParameters() when matrix has condensed slots.
* FIXED: Some functions now respect locally set mxOptions, when previously they would ignore them.

# OpenMx 2.6.9
* UPDATED: Stan header compatibility.

# OpenMx 2.6.8
* IMPROVED: Under the hood changes.

# OpenMx 2.6.7
* CHANGED: Default number of threads = 2. Previously OpenMx used (number of cores - 1)
	* This was done to reduce test-server burden for CRAN.
	* You can set threads using mxOption(NULL, "Number of Threads", cores)
	* nb: multi-threading supported on Linux and the OpenMx Team's build for Mac OS X.
* NEW: SLSQP multi-threading to evaluate the numerical gradient.
	* Use mxFitFunctionML(..., rowwiseParallel = FALSE).
	* By default with raw data mxFitFunctionML() parallelizes evaluation of the row likelihoods not the gradients.
* NEW: mxOption, "Parallel diagnostics". Set to "Yes", OpenMx provides diagnostic messages about the use of multiple threads.
* NEW: Functions dnbinom(), pnbinom(), dpois(), and ppois() (from the stats package) are now usable in MxAlgebras.
* IMPROVED: CSOLNP optimizer better at calculating CIs.
* IMPROVED: It is now possible to augment the GREML fit-function with an arbitrary scalar-valued function, to be evaluated and added to the fit function value. This can be used to regularize model-fitting with a prior log likelihood.
* IMPROVED: GREML fit function can also use analytic derivatives of the augmentation function.
* IMPROVED: mxRestore() now behaves correctly with argument strict = TRUE.
* FIXED: A subtle bug in the GREML fit function has been repaired.
	* Under certain circumstances, this bug caused analytic derivatives of the covariance matrix to fail to be recalculated after changes in the values of the free parameters upon which they depend.

# OpenMx 2.5.2
* NEW: mxFactorScores() enables regression factor-score estimates for RAM models!
* NEW: mxGenerateData() can generate data conditional on definition variables.
* NEW: SLSQP can use an analytic gradient during optimization.
* IMPROVED: mxTryHard() specially-tuned wrapper functions
	1. mxTryHardOrig()
	2. mxTryHardSSCT()
	3. mxTryHardWideSearch()
	4. mxTryHardOrdinal()
* NEW: imxRobustSE() calculates robust standard errors for parameter estimates from the sandwich estimator.
* NEW: functions in MxAlgebras:
	1. inverse trigonometric functions
	2. inverse hyperbolic functions
	3. logp2z() (standard-normal quantile function from log probabilities)
	4. lgamma1p() (accurate lgamma(x+1) for small x)
	5. Bessel functions, dbeta() pbeta()
* FIXED: Two mxGREMLDataHandler() behavior when blockByPheno = FALSE.
* FIXED: mxGREML automated handling of missing data when derivatives of the 'V' matrix are MxAlgebras.
* FIXED: LISREL path models now handle means correctly.
* FIXED: Factor-score estimates on factors with nonzero means resolved.
* IMPROVED: mxFactorScores() factor scores returned in original order
* IMPROVED mxFactorScores() no longer fails when SE are not available in the input model with type="ML" or type="WeightedML"
* IMPROVED: Several help pages have been updated, clarified, and made more complete.
* IMPROVED: Internal interface with NPSOL ensures optimizer consistently respects the "Major iterations" option.
* IMPROVED: Newton-Raphson optimizer handles encounters with parameter bounds less likely to cause convergence failure
* IMPROVED: mxGenerateData() now works with continuous-time state-space models.
* IMPROVED: Sufficient statistic likelihood adjusted to match the FIML value.
	* Prior versions (and Mx), did not correspond exactly to the full information formula.

# OpenMx 2.3.1
* NEW: Multi-group WLS has been implemented.
* NEW: Warning if Hessian is not convex at the solution (status code 5)
* NEW: mxRun() now displays the number of free parameters in the MxModel
* NEW: mxFactorScores() is now compatible with RAM models and multi-group models.
* NEW: coef() is now defined for MxModels (wrapper to omxGetParameters).
* NEW: mxCheckIdentification() is now compatible with GREML expectation.
* SPEEDUP: mxMatrix construction and modification
* SPEEDUP: NPSOL "warm start" now works correctly.
	* mxComputeGradientDescent() can provide optimizer with Cholesky of Hessian matrix at the start values. can reduce function evaluations.
* IMPROVED: mxTryHard().
* IMPROVED: mxGetExpected()'s compatibility with LISREL models
* IMPROVED: SLSQP ability when starting at a solution.
* IMPROVED: GREML matrix operations streamlined.
* IMPROVED: Evaluation of GREML analytic derivatives parallelized
	* Note this increases memory demand
* FIXED: A few bugs relating to state-space models have been repaired.
* FIXED: Serious GREML bugs fixed. Now safe to use mxFitFunctionGREML(dv=)
* NOTE: mxFactorScores() using 'ML' or'WeightedML' are deviations from the latent variable's mean.
	* If latent mean != 0, the scores must be shifted manually.


# Older changes

trunk
=====

* Prevent Varadhan2008 from failing near convergence
* Throw error on attempt to invert incomplete Hessian
* CSOLNP: Correct exclusion of inequality constraints from gradient/Hessian
* When checkpointing fit, record who requested it
* Restore parallel processing for CIs
* Add "Checkpoint Fullpath" to override output destination
* Make mxOption(model, val) return the global setting or the model's override
* Add WLS standard error and chi-square functions
* Consolidate starting value nudging logic in ComputeGD
* Try harder to keep CIs ordered properly vs estimate
* Add partial identification information to ID check.  Tells which
  parameters are not identified.
* CSOLNP: Avoid reuse of obm for more than 1 purpose
* CSOLNP: Simplify CI calculation
* CSOLNP: Tidy optimizer reporting
* CSOLNP: Count the number of major iterations
* Add Makefile rule to run the failing tests
* CSOLNP: Return stuff via CSOLNPFit; isolate Matrix inside of subnp
* CSOLNP: Remove some redundant information
* Remove template argument from CSOLNP::solnp
* Store equality & inequality results in CSOLNPFit
* Reduce variable lifetime
* Tidy some minor valgrind issues
* Remove q() from test
* CSOLNP: Convert control param to Eigen vector
* CSOLNP: Move bounds to CSOLNPFit; remove lots of deadcode
* CSOLNP: Move inequality bounds to CSOLNPFit
* CSOLNP: Remove Matrix arg from CSOLNPFit::solFun
* CSOLNP: Remove GLOB_ prefix since variables are no longer global
* CSOLNP: Eliminate more globals
* CSOLNP: Continue refactoring API
* CSOLNP: Pass parameter vector using Eigen
* CSOLNP: Rewrite omxProcessConstraintsCsolnp in Eigen
* CSOLNP: Don't bother creating the unused solEqB matrix
* CSOLNP: Eliminate reliance on all-zero solEqB
* Simplify setup of constraints
* CSOLNP: Put in temporary fix for conformability of derivs with constraints
* CSOLNP: Fix some problems revealed by Eigen conformability checking
* Rework calculation of compiler args; add IMX_SAFE switch
* Add test for omxGetParameters
* omxGetParameters: Omit duplicated parameters even when there are no submodels
* Tidy up MxMatrix-related error and warning messages; re-write MatrixErrorDetection.R
* BA81: Prohibit all NA rows
* Small but important MxMatrix-related adjustments.
* Add Jacobian-based model identification check.
* Added method for getting cov, means, and thresh from LISREL expectation.
* Added method for getting cov, means, and thresh from RAM expectation.
* Added generic for getting expected covariance, means, and thresholds.
* Small test for LISREL type, to be expanded
* LISREL model initialization.
* Better MxMatrix verification
* Conditionally "condense" the 'labels', 'free', 'lbound', and 'ubound' slots of MxMatrix objects.
* CSOLNP: Convert ind to Eigen; make easier to understand
* Enable ubsan for gcc 4.9+
* Make gcc happier
* CSOLNP: Remove more deadcode
* Avoid reuse of condif[123]
* test if OpenMx disorders thresholds
* CSOLNP: Remove bounds check; all Eigen ops are already bounds checked
* CSOLNP: Remove superfluous assignment
* CSOLNP: Remove deadcode
* CSOLNP: Converting function add to Eigen
* CSOLNP: Converting functions subtract and timess to Eigen
* CSOLNP: Converting functions subset, copy, multiply, and divide to Eigen
* CSOLNP: reapplying revisions 3976-3996
* Add gcc 4.9 support (same binary as gcc 4.8)
* WLS maybe working for joint.
* Helpful error when SaturatedLikelihood is given a list of two.
* Change the fake definition variable value for easier indexing
* CSOLNP: global variables moved to a struct
* Fix saturated DoF for raw data when only some of the variables are used.
* Bug fix for refs models with cov data
* CSOLNP: reverted version checked in
* Refrain from evaluating MxMatrix in the front-end
* Refactor computation of CFI, TLI, & RMSEA
* Move demos to nightly
* Prevent SEGV on ordinal columns with no thresholds
* Ensure continuous all-NA column is not run through mxFactor
* Add comment about multinomial degrees of freedom
* Enable --force-biarch for Windows fat binaries & remove obsolete rules
* various doc cleanup

Release 2.0.0-4004 (Oct 24, 2014)
=======================================

* Fix calculation of the saturated -2LL for IFA models
* Add some warning about CSOLNP
* Add hint about infeasible starting values (per Mike Neale)
* Rework capture of errors & warnings
* Permit running test from top dir
* CSOLNP: Unravel some confused Matrix resizing
* CSOLNP: Allow 0 coeff matrices; init to NaN instead of 1.0
* CSOLNP: Avoid allocation to report gradient
* CSOLNP: Add a variant of subset that copies out to Eigen::MatrixBase
* CSOLNP: Eliminate 2 instances of fillMatrix
* CSOLNP: Shorten lifetime of index variable
* CSOLNP: Convert a little bit of subnp to Eigen
* CSOLNP: For set{Row|Column}InPlace, treat 2nd arg as a vector
* Don't abort build if NPSOL is not found for a particular compiler
* CSOLNP: Switch Bcolj to stack allocation
* CSOLNP: Remove negate deadcode
* CSOLNP: Re-express findMax(matrixAbs()) as matrixMaxAbs()
* CSOLNP: Simplify divideByScalar2D & multiplyByScalar2D
* CSOLNP: Covert a few more copyInto -> copyIntoInplace
* CSOLNP: Move sob to stack allocated storage
* CSOLNP: Move alp to stack allocated storage
* CSOLNP: Move subnp_ctrl from Matrix to stack allocated Eigen::Array
* CSOLNP: Change setColumn to setColumnInplace
* CSOLNP: Change some setRow to setRowInplace
* Once constraints are turned into algebras, eliminate objectives from
  those new algebras and replace them with fitfunctions.
* CSOLNP: Convert some copyInto into copyIntoInplace
* CSOLNP: options added
* Update citation
* Drop psych dependency
* Document that mxLog will fail in Rgui on Windows
* Stop with error message if mxLog fails
* Updates to documentation
* Report CIs with NA as NA instead of throwing an error
* CSOLNP: void functions and memory allocation checks added to minor iteration loop
* Correctly collapse levels regardless of order (mxFactor)
* Repair derivatives of the beta density
* Avoid exception when only 1 optimizer has errors in the test suite
* Add mxFactor(..., collapse=TRUE)
* Fix bug that failed to print optimizer messages again in summary.  Add test.
* Fix mxEval(..., defvar.row) for when data is sorted
* Check for duplicated level names in factors in MxData
* Simple category collapsing for mxFactor
* Change mxStandardizeRAMpaths() to say if the model hasn't been run yet; amend test model accordingly.
* Add mxMakeNames & test
* Add eval by name wrapper for mxEval created by Spiegel
* 'make test' should record errors optimizer-wise
* Don't INSTALLMAKEFLAGS=""; allow setting from the environment
* Report 'make test' errors as they happen
* Throw warning from ref models when definition variables are present.
* Cope with zero length components in mxSimplify2Array
* Make mxFactor preserve rownames
* Update dims slot of expectations for models in mxRun.  Adjust saturated model accordingly.
* Make mxSimplify2Array preserve colnames
* Added threshold deviation labels for saturated models
* Fix bug in saturated model for raw matrix data.
* Implement qnorm() and lgamma() in mxAlgebras.
* Make models with row fitfunctions able to be re-run from an initial fit
  and end name collisions with filteredDataRow, existenceVector, and
  rowResults.
* Modify quoted formula error catch to only use match.call()$expression once.
* Smarter way to catch mxAlgebra formulas that are character strings.
* Catch quoted formula error in mxAlgebra
* Fix mxEval error and add test for square bracket with blank entry
* Modify wall time printing as per Dev mailing list discussion
* Small change so that mxStandardizeRAMpaths() handles independent submodels as intended.
* Adding correctly named Holzinger 1939 data set
* mxFitFunctionMultigroup: Fix error for unmatched submodel; add test
* omxQuotes should do something sensible when passed nothing
* Permit models with no matrices (to avoid masking other errors)
* Coerce x argument of mxFactor to character type
* Fix array indexing of ordinal thresholds (when there are unused columns in observed data)
* Add warning for data.type='cor'
* removing link to non-existant function in OpenMx [doc]
* Avoid stringification of non-finite numbers
* Charlie Driver pointed out that mxTryHard() fit attempts that result
  in npsolstatus -1 should be treated similarly to fit attempts ending
  in error
* Fix bug in summary for AIC/BIC with cor data.
* Drop duplicates from CI list
* Prevent acceleration of worsening fit
* BA81: Make E-step optimization less fragile

Release 2.0.beta3-3838 (Sep 26, 2014)
=======================================
* Constrain GRM item parameters more sensibly
* Newton-Raphson: Restore feasible parameter vector when fit obtains NaN
* Hopefully, fix performance decrement in FIML.
* Begin EM acceleration 1 cycle earlier
* Switch default EM acceleration to varadhan2008
* Reorg EM acceleration; add Varadhan & Roland (2008)
* Return SEs as "not requested" strings when argument SE=FALSE in mxStandardizeRAMpaths()
* Cope with observed data matrix with no NAs (observed data.frame not affected)
* Mention formulation of IFA independence model [doc]
* Drop support for non-double matrices
* CSOLNP: memory issue and reporting starting values as optimal parameters
* (Re?)enable parallel processing for confidence intervals
* Grab Eigen from RcppEigen
* Add ordinal ML/WLS test with both model identifications
* For saturated type=cov model, set TLI=1 and RMSEA=0
* Recompute thresholds (when provided)
* Rewrite algebra/fitfunction lookup in MxFitFunctionMultigroup
* Downgrade check of algebra dimnames to a warning for backward compatibility.
* Edit 'make clean' to also remove src/*.dll (Windows shared library).
* Add mxTryHard() and its man page.
* Add check for data type of algebra dimnames
* Implemented smarter tab-completion for MxObjects.
* Make the global freeGroup vector private
* Update SE study for relative tolerance
* SE simulation studies (IFA)
* Depend on package parallel (in R core since 2.14)
* Check threshnames for duplicates
* Store omxThresholdColumn in a std::vector
* LAD--CheckCode6.R moved to models/passing
* CSOLNP: status code report corrected
* Change the 'running model' output from mxRun to a message
* Toggle silent= back to FALSE as preference for running Reference/Saturated models.
* Add multilevel model example in state space form.  Data is grabbed from the web.
* Report optimizer in summary for default compute plan
* Rename nullModels -> refModels per dev discussion
* Add column of raw SEs to output of mxStandardizeRAMpaths(); edit its test model and man page appropriately.
* Fix 'unknown macro '\t'' warning I was getting from mxRestore man page
* Document Number of Threads option.
* Modify .svnignore files.  WLS branch and trunk changes.  Added
  continuous only WLS test.  Found and added state space algebra test.
* Include numObs in data for mxNullModels
* First draft of mxDataWLS
* Fixed documentation for reading checkpoint file via read.table
* Permit BA81 mean & variance specified with algebras
* Rename omxSaturatedModel to mxNullModels per dev discussion
* Refuse to clear model slots by assignment to NULL and say why
* Remove minItemsPerScore option
* Fix state space in FIML joint.

Release 2.0.beta2-3751 (Aug 20, 2014)
=======================================
* Switch default optimizer back to CSOLNP
* Made MxSummary have prettier printing of Chi-Square and RMSEA with CI.  Also added more checks for SummaryCheck.R
* Extinguish globalState.
* NLOPT and Simulated annealing added.
* Report condition number when standard errors are enabled.
* Improve destruction order
* Merge OMX_VERBOSE to OMX_DEBUG
* Simplify copyParamToModel API
* Move childList from omxState to FitContext
* Track whether a model has been run for submodels
* Modify FIML Single Iteration Joint to accommodate State Space expectations for continuous vars.
* Distinguish between stale models and models which has not been run
* Improve reporting of CIs for non-free parameters
* Ignore request for CIs if the label exists and is not a free parameter
* Add comment about RMSEA formula
* Reorganize omxSaturatedModel for better modularity
* Factor out processing of the run argument to omxSaturatedModel
* Remove our own version of std::max
* Consolidate FIML to omxFIMLSingleIterationJoint (except for state
  space). NOTE: This introduced a serious performance regression.
* Remove a bunch of duplicated code in FIML
* Add Oakes1999 EM information matrix method
* Set up setVarGroup handler for MxAlgebraFitFunction
* Add independent flag to mxComputeSequence
* Implement new unprotect strategy
* For summary(verbose=F), suppress SE and bounds if all NA
* minor error message change: show user the verboten list of names
* Various optimizations to improve frontend performance
* Don't blindly cast all NA columns to double. When the column is a
  factor, it needs to stay integer
* Only warn once about the model being modified since run
* Prohibit clashes between model and matrix names in type='RAM' models
* Clear modified-since-run flag on submodels
* Permit more than 1 fitfunction in mxComputeOnce
* The backend no longer resizes matrices that will be exported to R
* Switch to more accurate matrix log/exp
* Adjust pretty-printing of MxData to match slot names
* Print model name in summary
* Change error message in saturated model helper when algebra fit function detected.
* Improvements to mxSummary: RMSEA conf intervals, Information table
  adjustments, guidance to help page, statement about missing fit
  indices.
* CSOLNP: out of bound starting values issue solved. Mode added
* Revision to saturated model helper that allows multigroup and saturated models
* ComputeHessianQuality should not SEGV when no Hessian is available
* Add fix and test for mxFitFunctionR Hessian
* modified mxVersion to include R and platform + optimizer
* Add a few unprotects to reduce stack usage of big algebra
* Manual update from Hermine H Maes
* Extensive rework of the conformability checking pass
* Provide more helpful error messages when an expectation fails to initialize
* CSOLNP: potential fix for CI and non-linear constraint issues
* Add checks for AIC, BIC
* Add a conformability checking step to the backend
* better error when definition variable is missing
* Added names() support to the base R object types.  Changed the
  return value of MxModel names(). [For those who enjoy
  tab-completion]
* Make NPSOL the default if available
* Added verbose=FALSE argument to mxSummary to change amount of information that is printed.
* Add matrix logarithm algebra operator
* Move more CI code from frontend to backend
* Better way to label and reference figures (like APA style)
* Visually decorate figures so they stand out from the rest of the text [manual, HTML version]
* Recompute fit at the end of mxComputeGradientDescent
* CRAN prohibits variable length arrays
* Better error message when a CI is not found. Let user know what to check and do
* Fix for a minor bug that dates back to version 1.3.2 (at
  least). Now, mxVersion() reports the version number of the OpenMx
  package loaded into R's workspace, not the version number of the
  OpenMx package installed in the first directory in .libPaths()
* Added check of many (but not all) of the parts of MxSummary.
* Added joint ordinal continuous example and an example with means to LISREL man page.
* Improve safety of mxFactor
* Automatically extract and run code from the manual to ensure that it works
* Added reporting of Chi square DoF.
* Minor bug fix in saturated model helper.
* Fixed bug in Chi square degrees of freedom calculation.
* modified error message in mxRun and mxOption to tell user how to change
  default optimizer, and how to add a fit function
* Refuse to set "Default optimizer" on models
* Reinstate R_CheckUserInterrupt (got commented out by mistake)
* Fail if a model has an expectation, no fitfunction, and no custom compute plan
* Fix signature mismatch in displayCompute
* Remove non-ascii em dashes
* Fix our signature for logLik S3 method
* If an MxMatrix is constant, avoid copy
* Advertise our version as 2.0.0 instead of 999.0.0
* Re-implement omxData storage to facilitate dynamic data
* omxAssignFirstParameters should not explode with 0 free parameters
* Move packageStartupMessage to .onAttach
* Reject vector=TRUE as part of MxFitFunctionMultigroup
* Make FIML respect vector=FALSE

Release 2.0.beta1-3473 (May 30, 2014)
=======================================
* Fix for R 3.1
* Add packageStartupMessage if compiled without OpenMP
* Remove extra copy of # of evaluations
* Fix algebra dependency tracking
* Name anonymous algebras to aid debugging
* numThreads is always 1 without OpenMP
* MxRAMModel should not assume MxExpectationRAM
* Permit optimization directly on mxFitFunctionML(vector=TRUE)
* Remove unimplemented and crashing omxSetFinalReturnsFIMLFitFunction
* Added a warning about using "one" as a label in mxPaths.
* minor error message mod to suggest action if non-square matrix declared as cov
* Warn if summary is called on a model that was modified after mxRun
* Bug Fix: the backend resized 1x1 matrices to MxN when used in scalar/elementwise multiplication, but then did not resize them back in the frontend.
* Fix mxRename with constraints bug.
* Learn mxMatrix dimnames from values if explicit dimnames omitted
* fix issue 'summary() of fitted mxModel object returns error'
* Confidence Intervals now give the name of the parameter if they are from an MxMatrix instead of the MatrixName[row, col].
* Exterminate strncmp
* Distinguish between sorted/unsorted and whether sorting is requested
* Remove unneeded parameter from isErrorRaised()
* Fixed Saturated Likelihood Bug by not populating that attribute for FIML.
* Support NPSOL warmStart
* Add option to checkpoint every evaluation
* Bug fix the number of observed statistics in non-IFA models with raw data.
* In addition to the mxOption, let mxData(sort=FALSE) request unsorted data
* Add a git bisect script
* Fixed bug where mxOption 'Default optimizer', and other options as well, get overwritten to the global defaults
* logLik for mxModel (contributed by Andreas Brandmaier)
* mxOption without a value to show the current setting
* Place likelihood-based CI code into a separate MxCompute step
* Add rownames to standardErrors
* Do not resetDataSortingFlags
* Remove unused arguments from omxRaiseError, omxRaiseErrorf
* Remove unused argument from omxResizeMatrix
* Handle omxRemoveRowsAndColumns for row-major order; add tests
* omxAliasMatrix, omxResetAliasedMatrix are essentially omxCopyMatrix
* Simplify implementation of aliased matrices
* Fix memory corruption in omxResetAliasedMatrix
* Checkpoint rewrite
* Allow standard errors without Richardson extrapolation
* Preserve dimnames in assignments to MxMatrix slots
* Rework reporting of iterations and optimizer status; add failing test
* Fix memory corruption in matrix populate lists
* Remove unused omxInitMatrix argument
* Add more gdb hints
* Fix misuse of stdargs macros
* Enable gdb for the regular 'make test' variants
* Hook up recordIterationError for FIML
* Add FitContext::recordIterationError
* Enumerate errors if there is more than one
* Report all the errors instead of only the most recent one
* mxCI should not interfere with SEs (or fit)
* Add dimnames to @output$hessian and @output$ihessian
* Fix omxMatrix memory model
* Remember how many cores we detected
* Fix mxCI for a vector of parameter names; add test
* Hook up numerical integration precision parameters
* Improve detection of CPU architecture
* Minor bug fix in saturated model generator
* Pull libnpsol.a from our website (linux only, so far)
* Hint the correct way to customize the compiler/compiler options
* Add info on using omxAssignFirstParameters with SetParameters
* Add mxOption for max stack depth
* Factor out lots of calls to Rf_mkChar
* Remove "at iteration" from error messages (simplifies log diff)
* Make elementwise algebra ops conformable for the scalar-matrix case
* Continue with omxMatrix API simplification
* Enable R_NO_REMAP for a cleaner namespace
* Add check for algebra dimnames
* Don't use R to calculate our algebra result matrix dimensions
* Initial compute protocol for the whole tree of dependencies
* Don't rely on R to evaluate our algebra
* Remove non-reproducible pointer addresses from logs (makes it possible to diff logs)
* Change backend initialization order
* Set up the usual gdb breakpoints automatically
* Update references to mxFitFunctionAlgebra and correct dot multiplication explanation in mxAlgebra.Rd
* Improve some conformability checks by showing the dimension mismatch
* Always report @output$fit but don't report @output$minimum unless it is
* Using the name "stderr" with #pragma omp critical can cause conflicts when building multi-thread binaries under Windows.
* Permit rescaling of log-likelihood
* Compute condition number of information matrix
* Control some ComputeEstimatedHessian knobs from R
* Do not use roxygen to create Collate field
* Reduce use of protect stack
* Suggest how to debug protect stack overflow
* Permit easy specification of the default optimizer in R code
* Keep parameters within bounds
* Rename most PPML functions to keep them un-exported
* Allow summary to work with unnamed estimates
* Fixed a bug where ML fit functions with raw data (FIML) and vector=TRUE set were returning a single value instead of a vector.
* Updating twinData.Rd to document the reuse of zyg 6:10. Also added Nick Martin reference
* Document how to find & adjust R's default compiler flags
* Add back some UNPROTECTs
* Be more paranoid about importing mxData to backend
* Avoid final copyParamToModel when not needed
* Remove obsolete performance counters
* Fix/remove improper printf style formats
* Change how dirty matrices are indicated
* Rewrite dependency tracking
* Grab expectation names in the backend
* Disable NPSOL gradient verification by default (this is a developer feature)
* Permit better control over whether NPSOL uses gradients
* Remove most instances of setFinalReturns, mark deprecated
* Don't communicate unused dependency information to the backend
* Store fit function name (regression fix?)
* Free parameter groups
* Store omxGlobal as a pointer to ensure proper init & destruction
* Threads don't help RFitFunction, force single threaded
* NPSOL linear constraints are unused (deadcode)
* Set up a build rule for installing without NPSOL
* Replace bloated autoconf script with a simple shell script
* Prevent developers from using Rprintf
* Split omxState into truly global stuff and per-thread stuff
* Modified mxSaturatedModel to take either a model or a data set, and added a run= argument with FALSE as the default.
* mxSaturated model support for ordinal data.
* Skip recomputing thresholds when they don't exist
* Replace completely broken PBS cpus detection
* Add thread-safe logging functions
* Changed all objective functions to return NAMED lists of expectation and fitfunction.
* Switch over to new structured Compute system
* Added manual page for mxFitFunctionML
* Encapsulate NPSOL into omxComputeGD; rework error status reporting
* Pass blank gradient and hessian to omxNPSOLConfidenceIntervals
* Reorder matrix name processing
* Dump matrices using valid R syntax (debugging)
* Remove unimplemented OMX_SOCKET_CHECKPOINT
* Simplify management of omxData.type
* Switched backend from C to C++
* Fixed maximum number of dimensions bug so that only dimensions other than those -inf to +inf are counted
* Add API for internal expectations; remove redundant copy of rObj
* C-side support for submodels
* Add deps on MASS and mvtnorm
* Remove lots of UNPROTECTs
* Refactor allocation of omxAlgebra.args
* Teach MxExpectation to store its relationship in the model tree
* Add default A,S,F matrix names for mxExpectationRAM
* Process submodels concurrently with manifest & latent variables (early)
* Simplify error reporting of unrecognized arguments to mxModel
* Factor out interpretation of mxPath's connect argument
* fix insertMeansPathRAM to catch arrows=2
* minor fix to insertMeansPathRAM() 1. catch arrows != 1 in means paths; 2. reword error
* Automatically balance PROTECT/UNPROTECT
* Added suggestion for solution when mxData are cov/cor and not symmetric.
* Improved mxData error messages
* Added cov2cor and chol to supported functions.
* Update mxCompare() to handle missing comparison parameter
* update showFitStatistics() to handle empty compareSummaries list
* Added comment to RAM model error, prompting user to view ?mxData when their model has no data
* Add varargs replacement for omxRaiseError
* Make omxCheckCloseEnough compare missingness pattern too
* Complain if any starting values are missing
* Fix uninitialised memory access in omxSelectRows & omxSelectCols,
  Bug only affected non-square matrices
* changed FIML fit functions to report all elements when vector=TRUE
* added mxThreshold() function
* migrated from MxObjective* to mxFitFunction* and mxExpectation*
* bug fix in multiple-iteration upper confidence interval estimation
* bug fix in mxRun() where a model with grandchildren submodels can
  have its internal state corrupted
* Added vech2full and vechs2full mxAlgebras: inverses of vech and vechs
* detect non-conformable arrays for elementwise division in the front-end
* fix crash in Hessian calculation when no objective function is specified
* warn if standard errors are enabled and the numeric Hessian calculation is disabled
* added confidence interval calculations to checkpoint output

Release 1.3.0-2168 (September 17, 2012)
=======================================
* added mxOption() for "Analytic Gradients" with possible values "Yes"/"No"
* added 'cache' and 'cacheBack' arguments to mxEval()
* added omxLocateParameters() function (see ?omxLocateParameters)
* type='RAM' allowing 'manifestVars' argument to appear in a different order
  than in the observed covariance matrix.
* the configuration mxRun(model, checkpoint = TRUE) writes a line
  in the checkpoint file at the conclusion of model optimization.
* added "Always Checkpoint" to mxOptions() with values "Yes" or "No"
* header for the checkpoint file will identify anonymous 
  free parameters with the string modelName.matrixName[row,col]
* omxGetParameters() and omxSetParameters() support anonymous
  free parameters
* implemented cov2cor in the OpenMx backend
* mxOption "Major iterations" accepts either a value or a function
* now tracking the MxAlgebra and MxMatrix objects that need to be updated
  when populating free parameters or updating definition variables.
* performance improvements in mxModel() when building RAM models
* performance improvements in mxRun() frontend for large matrices
* rewrite on the processing of objective functions in the frontend
* added R functions omxCbind(), omxRbind(), and omxTranspose()
* added 'fetch' argument to omxGetParameters()

Release 1.2.5-2156 (September 5, 2012)
======================================
* bugfix to omxRAMtoML() when input model has covariance data
* fixing memory leak in cleaning up algebras when optimization is complete
* fixing memory leak in omxImaginaryEigenvalues()
* fixed a bug that manifests in confidence intervals with definition variables
  (see http://openmx.ssri.psu.edu/thread/1505)
* fixed a bug in the identification of NA definition variables
  (see http://openmx.ssri.psu.edu/thread/1521)

Release 1.2.4-2063 (May 22, 2012)
=================================
* added argument "name" to omxSetParameters()
* error checking for 0-length arguments to mxPath()
* fixing several memory leaks

Release 1.2.3-2011 (April 10, 2012)
===================================
* bugfix when (I-A)^-1 speedup is disabled
* bugfix in mxAlgebra() detection of missing operator or function
* bugfix in joint FIML when: 1) there are no definition variables, 
  2) the first (sorted) row with a new number of ordinal variables, and 
  3) has all continuous variables missing

Release 1.2.2-1986 (March 22, 2012)
===================================
* setting default value for mxOption("UsePPML") to "No" until
  the feature is adequately tested. Added test enforcing default
  value to test suite.

Release 1.2.1-1979 (March 21, 2012)
===================================
* bug fix for interaction of matrix transpose and square bracket substitutions
* renaming mxLISRELObjective() to imxLISRELObjective(), LISREL objective function is not yet implemented.
* improved error messages when a free parameter has multiple lbounds/ubounds
* improved error messages when passing strings into mxModel()
* improved error messages in mxEval()

Release 1.2.0-1926 (February 04, 2012)
======================================
* new interface such that 'objective@info$expMean', 'objective@info$expCov'
  and 'objective@info$likelihoods' are available in FIML or RAM objective functions.
* bug fix for omxSetParameters with lbound or ubound of NA.
* catching unknown matrix operator or function in mxAlgebra() or mxConstraint() declaration
* deprecated 'all' argument from mxPath function and replaced it with 'connect' argument. Also updated demos that used 'all' argument and documentation.
* removed 'excludeself' argument from mxPath() function
* added deprecation warning for argument 'all' = TRUE in mxPath()
* added support for OpenMP in hessian calculation
* added support for OpenMP in continuous FIML objective function
* added support for OpenMP in ordinal FIML objective function
* added support for OpenMP in joint ordinal/continuous FIML objective function
* added more descriptive error message to mxOption()
* serializing the sum operation on likelihoods - floating point addition is not associative

Release  1.1.2-1818 (October 24, 2011)
======================================
* fixed bug in summary() function, see http://openmx.ssri.psu.edu/thread/1104
* fixed bug in independence model likelihood calculation for ML objective
* included support for gcc 4.5 and 4.6 under 32-bit x86

Release  1.1.1-1784 (September 11, 2011)
========================================
* fixed several bugs in joint ordinal-continuous integration
* fixed several typos in User Guide

Release  1.1.0-1764 (August 22, 2011)
=====================================
* Joint and Ordinal documentation included and up to date.
* added deprecation warning for argument 'all'=TRUE in mxPath()
* added omxSelectRowsAndCols, omxSelectRows and omxSelectCols documentation
* fixed model flattening to keep track of confidence intervals in submodels
* added error checking for dimnames on MxMatrix and MxAlgebra objects
* reformatted comments and heading style for all demos
* fixed segmentation fault on backend error condition
* turn off jiggling of free parameters with starting values of 0.0 when useOptimizer=FALSE
* allow non-RAM objective functions in RAM model
* change model type names to 'default' and 'RAM'
* no longer explicitly transforming RAM + raw data models into FIML models
* fixed bug in MxMatrix indexing operator
* added argument 'threshnames' to mxFIMLObjective() and mxRAMObjective()
* renaming majority of omx* functions to imx* functions. See http://openmx.ssri.psu.edu/thread/761
* added error checking to mxOption() function
* added "Optimality tolerance" to mxOption() selection
* added 'lbound' and 'ubound' columns to summary() output of free parameters
* added asterisks to the 'lbound' and 'ubound' columns when feasibility tolerance is met
* added "omxNot" function to the set of available mxAlgebra() function
* added "omxSelectRows", "omxSelectCols", and "omxSelectRowsAndCols" as mxAlgebra operators()
* added "mean" function to the set of available mxAlgebra() functions
* added slots "expCov" and "expMean" to the MxRAMObjective function
* added useOptimizer option to mxRun.
* added error checking in frontend and backend for non-positive-definite observed covariance matrices
* added "omxGreaterThan", "omxLessThan", "omxApproxEquals", "omxAnd", and "omxOr" operators to the set of mxAlgebra() operators
* error checking for model[[1]] or model[[TRUE]]
* error checking in the front-end whether more than 20 ordinal columns are present in a data set
* improved performance in the front-end in mxModel() for adding paths to RAM models
* print name of algebra when operator has too few or too many arguments
* added mxErrorPool() function and R documentation.
* added Apache license information to all R documentation files.
* new implementation of mxEval().
* new argument 'defvar.row' to mxEval().  See ?mxEval.
* handling definition variables for (I - A) ^ - 1 speedup
* handling square bracket labels for (I - A) ^ - 1 speedup
* added argument 'free' to omxGetParameters.  See ?omxGetParameters.
* added argument 'strict' to omxSetParameters.  See ?omxSetParameters.
* eliminated warnings for confidence interval optimization codes
* added "..." argument to mxRObjectiveFunction()
* fix memory leak in RAM objective function
* removed dependency to MBESS library in R documentation
* added more descriptive error message when thresholds are not sorted
* incorporated NaN unsafe matrix-matrix multiplication (dgemm) from R <= 2.11.1
* incorporated NaN unsafe matrix-vector multiplication (dgemv) from R <= 2.11.1
* return NA in mxVersion() if "OpenMx" cannot be found
* fix infinite loop in objective function transformations
* added initialization to load OpenMx on swift workers
* implemented omxParallelCI() to calculate confidence intervals in parallel
* only calculating CIs for upper triangle of symmetric matrices
* cleanup appearance of transient MxMatrix objects in error messages
* fix bug with very large number of omxUntitledName() objects
* added optional argument CPUS=n to "make test" target
* change snowfall interface to use sfClusterApplyLB()
* storing raw data in row-major order, and copying contiguous data rows
* fix bug in mxModel() when using remove = TRUE

Release 1.0.7-1706 (July 6, 2011)
=================================
* error checking in front-end for non-positive-definite observed covariance matrices
* fix bug in MxMatrix indexing operator
* added deprecation warning for argument 'all'=TRUE in mxPath()

Release 1.0.6-1581 (March 10, 2011)
===================================
* Bug fix corner case in sorting data with definition variables

Release 1.0.5-1575 (March 8, 2011)
==================================
* added error checking in ML objective to match expected covariance and observed covariance matrices
* updated BootstrapParallel.R demo to use mxData() instead of model@data
* added dataset from the Psychometrika article <www.springerlink.com/content/dg37445107026711>

Release 1.0.4-1540 (January 16, 2011)
=====================================
* added initialization to load OpenMx on swift workers
* only calculating CIs for upper triangle on symmetric matrices
* fix bug with very large number of omxUntitledName() objects
* incorporated NaN unsafe matrix-vector multiplication (dgemv) from R <= 2.11.1
* fix bug in mxModel() when using remove = TRUE
* Added intervals to MxModel class documentation

Release 1.0.3-1505 (November 10, 2010)
======================================
* return NA in mxVersion() if "OpenMx" cannot be found
* eliminate infinite loop in objective function transformations

Release 1.0.2-1497 (November 5, 2010)
=====================================
* missing <errno.h> include 
* fix memory leak in RAM objective function
* incorporated NaN unsafe matrix-matrix multiplication (dgemm) from R <= 2.11.1

Release 1.0.1-1464 (October 8, 2010)
====================================
* bugfix for mxEval() and MxData objects
* handling definition variables for (I - A) ^ - 1 speedup
* handling square bracket labels for (I - A) ^ - 1 speedup
* added argument 'free' to omxGetParameters.  See ?omxGetParameters.
* added argument 'strict' to omxSetParameters.  See ?omxSetParameters.
* eliminated warnings for confidence interval optimization codes

Release 1.0.0-1448 (September 30, 2010)
=======================================
* added missing entries to demo 00INDEX file

Release 0.9.2-1446 (September 26, 2010)
=======================================
* added growth mixture models to user guide
* added initial Swift hook in omxLapply() - currently activated only for mxRun calls
* fixed a bug in mxRename() when encountering symbol of missingness
* bugfix for crash when '*' is used instead of '%*%'
* feature removal: square brackets in MxMatrix labels now accept only literal values
* bugfix for definition variables used in mxRowObjective (which is still experimental)
* bugfix for (I - A) ^ -1 speedup with FIML optimization

Release 0.9.1-1421 (September 12, 2010)
=======================================
* fixed a bug in -2 LL calculation in RAM models with definition variables and raw data

Release 0.9.0-1417 (September 10, 2010)
=======================================
* improved error messages for non 1 x n means vectors in FIML and ML
* fixed a performance bug that was forcing too many recalculations of the covariance matrix in FIML optimizations.
* default behavior is to disable standard error calculations when model contains nonlinear constraints
* improved error messages for NA values in definition variables
* added error message when expected covariance dimnames and threshold dimnames do not contain the same elements.
* fixed a bug when mxRename() encounters a numeric or character literal.
* new chapters added to OpenMx user guide.
* fixed a bug in -2 log likelihood calculation with missing data

Release 0.5.2-1376 (August 29, 2010)
====================================
* improved error messages when identical label is applied to free and fixed parameters
* added 'onlyFrontend' optional argument to mxRun() function.  See ?mxRun.
* disabling cbind() and rbind() transformations as they are broken.

Release 0.5.1-1366 (August 22, 2010)
====================================
* added error detection when multiple names are specified in mxMatrix(), mxAlgebra(), etc.
* removed 'digits' argument from mxCompare.  Target behavior of argument was unclear.  See ?mxCompare
* more informative error messages for 'dimnames' argument of objective functions
* more informative error messages when constraints have wrong dimensions
* improved error detected for 'nrow' and 'ncol' arguments of mxMatrix() function
* fixed a bug in ordinal FIML objective functions with non-used continuous data

Release 0.5.0-1353 (August 08, 2010)
====================================
* calculating cycle length of RAM objective functions
* bugfix: preserve rownames when converting data.frame columns to numeric values
* 'nrow' and 'ncol' arguments now supersede matrix dimensions in mxMatrix()
* add boolean argument 'vector' to mxRAMObjective() for returning the vector of likelihoods
* added demo(OneFactorModel_LikelihoodVector) as example of 'vector=TRUE' in RAM model
* cbind() and rbind() inside MxAlgebra expressions with all arguments as MxMatrix objects are themselves transformed into MxMatrix objects
* bugfix with square bracket substitution
* finishing implementing sorting of raw data in mxFIMLObjective()
* added 'RAM Optimization' and 'RAM Max Depth' to model options.  See ?mxOption
* added support for linux x86_64 with gcc 4.1.x

Release 0.4.1-1320 (June 12, 2010)
==================================
* confidence interval optimizations now jitter if they can't get started
* added error checking for dimensions of expected means in ML + FIML objectives
* check for missing observed means when using (optional) expected means in ML
* added citation("OpenMx") information

Release 0.4.0-1313 (June 09, 2010)
==================================
* fixed bug in calculation of confidence intervals around non-objective values
* checking for partial square bracket references on input
* fixed error reporting for non-positive-definite covariances in FIML
* implemented initial sorting-based speedup for FIML objectives
* added 'No Sort Data' to mxOptions()
* eliminated getOption('mxOptimizerOptions') and getOption('mxCheckpointOptions') 
* added getOption('mxOptions')
* error messages for illegal names provide function call information
* added 'estimates' column to confidence intervals in summary() of model
* fixed bug so checkpointing will work in R 2.9.x series
* fixed bug so mxRename() works on confidence interval specification
* renaming 'estimates' column of confidence interval summary output to 'estimate'

Release 0.3.3-1264 (May 24, 2010)
=================================
* confidence interval frontend was requesting nonexistent matrices
* omxNewMatrixFromMxMatrix() assumed input was always integer vector S-expression
* confidence intervals mislabeling free parameter names

Release 0.3.2-1263 (May 22, 2010)
=================================
* added 'newlabels' argument to omxSetParameters() function
* now throwing errors to the user when detected from the backend in mxRun()
* checkpointing mechanism implemented - mxRun(model, checkpoint = TRUE)
* never computing confidence intervals for matrix cells where free = FALSE (all bets are off on algebras)
* added mxOption(model, "CI Max Iterations", value) 
* added documentation for mxRestore() function
* default expected means vectors are no longer generated

Release 0.3.1-1246 (May 09, 2010)
=================================

* new arguments to mxRun() for checkpointing and socket communication (doesn't work yet)
* throw an error if FIML objective has thresholds but observed data is not a data.frame object.
* bugfix for R version 2.11.0 is detecting as.symbol("") character as missing function parameter
* mxFactor() function accepts data.frame objects
* added mxCI() function to calculate likelihood-based confidence intervals
* mxCompare() shows model information for base models
* added flag all=[TRUE|FALSE] to mxCompare() function
* 'SaturatedLikelihood' argument to summary() function will accept MxModel object

Release 0.3.0-1217 (Apr 20, 2010)
=================================

new features
------------
* implemented new ordinal data interface http://openmx.ssri.psu.edu/thread/416#comment-1421 (except for 'means=0' component)
* added mxFactor() function (see ?mxFactor for help)
* added R documentation for rvectorize() and cvectorize()
* implemented eigenvalues and eigenvectors
* added 'numObs', 'numStats' arguments to mxAlgebraObjective()
* added 'numObs', 'numStats' arguments to summary()

changes to interface
--------------------
* mxConstraint("A", "=", "B") is now written as mxConstraint(A == B)
* renaming 'cov' argument to 'covariance' in omxMnor()
* renaming 'lbounds' argument to 'lbound' in omxMnor()
* renaming 'ubounds' argument to 'ubound' in omxMnor()
* renaming 'cov' argument to 'covariance' in omxAllInt()
* argument 'silent = TRUE' to mxRun() will no longer suppress warnings
* argument 'suppressWarnings = TRUE' to mxRun() will suppress warnings

bug fixes
---------
* added error checking for mxBounds() with undefined parameter names
* improving error messages in mxMatrix() function
* fixed aliasing bug in mxAlgebra() with external variables and constant values

internal
--------
* added directory to repository for nightly tests ("make nightly")
* performance improvement to namespace conversion
* changing MxPath data structure from a list to an S4 object
* new signature for omxSetParameters(model, labels, free, values, lbound, ubound, indep)
* checked in implementation of merge sort
* changed mxCompare() signature to mxCompare(base, comparison, digits = 3)


Release 0.2.10-1172 (Mar 14, 2010)
==================================
* bugfix for assigning default data name when default data does not exist
* bugfix for sharing/un-sharing data to a three-level hierarchy
* bugfix for error reporting in bad matrix access, un-calculated std. errors, and poor omxAllint thresholds
* implemented rvectorize, cvectorize algebra functions: vectorize by row, and vectorize by column
* not allowing the following forbidden characters in names or labels: "+-!~?:\*/^%<>=&|$"
* added timestamp to summary() output
* re-implemented summary() function to handle unused data rows, and independent submodels
* re-implemented names(model), model$foo, and model$foo <- value to return all components of a model tree
* started work on an "OpenMx style guide" section to the User Guide.
* fix documentation + demo errors brought to our attention by @dbishop

Release 0.2.9-1147 (Mar 04, 2010)
=================================
* bugfix for ordinal FIML with columns that are not threshold columns
* bugfix for detection of algebraic cycles when multiple objective functions are present
* test cases included for standard error calculation
* enabled standard error calculation by default

Release 0.2.8-1133 (Mar 02, 2010)
=================================
* bugfix for memory leak behavior in Kronecker product calculation
* implemented Kronecker exponentiation operator %^%.
* actually updating means calculations in ordinal FIML models
* removing standard errors from summary() until they are computed correctly

Release 0.2.7-1125 (Feb 28, 2010)
=================================
* added 'unsafe' argument to mxRun function.  See ?mxRun for more information.
* bugfix for filtering definition variable assignment to current data source.
* bugfix for using data.frame with integer type columns that are not factors.

Release 0.2.6-1114 (Feb 23, 2010)
=================================
* implemented omxAllInt, use ?omxAllInt for R help on this function
* implemented omxMnor, use ?omxMnor for R help on this function
* added option to calculate Hessian after optimization.  See ?mxOption for R help.
* added option to calculate standard errors after optimization.  See ?mxOption for R help.
* added R documentation for vech, vechs, vec2diag, and diag2vec
* performance improvements for independent submodels
* added 'silent=FALSE' argument to mxRun() function
* added R documentation for omxApply(), omxSapply(), and omxLapply()
* wrote mxRename() and added R documentation for function
* added 'indep' argument to summary() to ignore independent submodels (see ?summary)
* added 'independentTime' to summary() output. Wall clock time for independent submodels.
* added 'wallTime' and 'cpuTime' to summary() output. Total wall clock time and total cpu time.
* implemented ':' operator for MxAlgebra expressions. 1:5 returns the vector [1,2,3,4,5]
* implemented sub-ranges for '[' operator in MxAlgebra expressions. foo[1:5,] is valid inside algebra.
* added omxGetParameters, omxSetParameters, omxAssignFirstParameters.  Use ? for documentation.
* renamed all objective function generic functions from omxObj* to genericObj*
* enumerated OpenMx and NPSOL options in ?mxOption documentation
* implemented foo[x,y] and foo[x,y] <- z for MxMatrix objects


Release 0.2.5-1050 (Jan 22, 2010)
=================================
* added 'mxVersion' slot to output of summary() function.
* added documentation of summary() function.
* set default function precision for ordinal FIML evaluation to "1e-9"
* not throwing an error on mxMatrix('Full', 3, 3, labels = c(NA, NA, NA))
* applying identical error checking to single thresholds matrix and single thresholds algebra
* implemented generic method names() for MxModel objects.
* implemented vec2diag() and diag2vec() matrix algebra functions.

Release 0.2.4-1038 (Jan 15, 2010)
=================================
* definition variables can now be used inside algebra expressions
* definition variables inside of MxMatrices will populate to the 1st row before conformability checking. In plain english: you do not need to specify the starting values for definition variables.
* the square-bracket operator when used in MxMatrix labels is no longer restricted to constants for the row and column.  The row and column arguments will accept any term that evaluates to a scalar value or a (1 x 1) matrix.
* summary() on a model returns a S3 object.  Behaves like summary() in stats package.
* eliminated UnusualLabels.R test case.  Too many problems with windows versus OS X versus linux.
* implemented vech() and vechs() functions: half-vectorization and strict half-vectorization
* fixed bug in ordinal FIML when # of data columns > # of thresholds
* added 'frontendTime' and 'backendTime' values to summary() output. They store the elapsed time of a model in the R front-end and C back-end, respectively.
* created a name space for the OpenMx library. Only mx**() and omx**() functions should be exported to the user, plus several miscellaneous matrix functions and S4 generic functions.
* corrected 'observedStatistics' output of summary() to exclude definition variables
* corrected 'observedStatistics' output of summary() count the number of equality constraints

Release 0.2.3-1006 (Dec 04, 2009)
=================================
* added 'vector' argument to mxFIMLObjective() function. Specifies whether to return the likelihood vector (if TRUE) or the sum of log likelihoods (if FALSE). Default value is FALSE.
* renamed omxCheckEquals() to omxCheckIdentical(). omxCheckIdentical() call "identical" so that NAs can be compared.
* added checking of column names of F and M matrices in RAM objective functions.
* added 'dimnames' argument to mxRAMObjective() function. Populates the column names of F and M matrices.
* added square bracket operator to MxAlgebra expressions.  A[x,y] or A[,y] or A[x,] or A[,] are valid.
* square bracket operator supports row and column string arguments.
* mxModel(remove=TRUE) accepts both character names or S4 named entities.
* added support for x86_64 on OS X 10.6 (snow leopard)
* fixed support for x86_64 on Ubuntu 9.10 (gcc 4.4)
* throw error message when inserting a named entity into a model with an identical name
* added Anthony William Fairbank Edwards "Likelihood" (1972; 1984) A, B, O blood group example to online documentation

Release 0.2.2-951 (Oct 29, 2009)
================================
* omxGraphviz() either prints to stdout or to a filename
* updated omxGraphviz() to draw an arrow if (value != 0 || free == TRUE || !is.na(label))
* omxGraphviz() returns a character string invisibly
* error checking for bogus definition variables
* summary() uses matrix dimnames by default, use options('mxShowDimnames'=FALSE) to disable
* added support for gcc 4.4 (Ubuntu 9.10) on x86 and x86_64 architectures
* created R documention for omxGraphviz() function
* generalized dependency specification for objective functions (omxObjDependencies)
* fixed cross-reference links in User Guide

Release 0.2.1-922 (Oct 10, 2009)
================================
* checked observed data for dimnames on ML objective functions
* using dimnames= argument to mxMLObjective() and mxFIMLObjective()
  propagates to MxMatrix objects on output.
* bug fix for error message where model name is incorrect
* updated user guide in response to feedback from beta testers
* incremented version number to 0.2.1-922 to sync demos and user guide

Release 0.2.0-905 (Oct 06, 2009)
================================
* several of the twin model demo examples have been re-written.
* fixed bug with non-floating point matrices.
* more error checking for mxPath().
* tools/mxAlgebraParser.py will convert Mx 1.0 algebra expressions (Python PLY library is required).
* renamed "parameter estimate" column to "Estimate" and "error estimate" to "Std.Error" in summary().
* added 'dimnames' argument to mxFIMLObjective() and mxMLObjective().
* error checking for RAM models with non-RAM objective functions.

Release 0.1.5-851 (Sep 25, 2009)
================================
* improved error messages on unknown identifier in a model (beta tester issue)
* fixed bug in mxMatrix() when values argument is matrix and byrow=TRUE
* implemented square-bracket substitution for MxMatrix labels
* fixed a bug in computation of omxFIMLObjective within an algebra when definition variables are used
* significant alterations to back-end debugging flags
* tweaked memory handling in back-end matrix copying
* added support for x86_64 linux with gcc 4.2 and 4.3

Release 0.1.4-827 (Sep 18, 2009)
================================
* added checking and type coercion to arguments of mxPath() function (a beta tester alerted us to this)
* moved matrices into submodels in UnivariateTwinAnalysis_MatrixRaw demo
* added Beginners Guide to online documentation
* mxRun() issues an error when the back-end reports a negative status code
* named entities and free or fixed parameter names cannot be numeric values
* constant literals are allowed inside mxAlgebra() statements, e.g. mxAlgebra(1 + 2 + 3)
* constant literals can be of the form 1.234E+56 or 1.234e+56.
* type checking added to mxMatrix arguments (prompted by a forum post)
* mxPath() issues an error if any of the arguments are longer than the number of paths to be generated
* data frames are now accepted at the back-end
* FIML ordinal objective function is now working. Still a bit slow and inelegant, but working
* FIML ordinal now accepts algebras and matrices. dimnames of columns must match data elements
* implemented free parameter and fixed parameter substitution in mxAlgebra statements
* implemented global variable substitution in mxAlgebra statements
* turned off matrix and algebra substitution until a new proposal is decided
* snow and snowfall are no longer required packages
* added cycle detection to algebra expressions
* mxEval() with compute = TRUE will assign dimnames to algebras
* added dimnames checking of algebras in the front-end before optimization is called
* added 'make rproftest' target to makefile

Release 0.1.3-776 (Aug 28, 2009)
================================
* mxEvaluate() was renamed to mxEval() after input from beta testers on the forums.
* new function mxVersion that prints out the current version number (beta tester request).
* When printing OpenMx objects, the @ sign is used where it is needed if you would want to print part of the object (beta tester request).
* now supports PPC macs.
* implemented AIC, BIC and RMSEA calculations.
* mxMatrix documentation now talks about lower triangular matrices (beta tester request).
* fixed bugs in a number of demo scripts.
* added chi-square and p-value patch from beta tester Michael Scharkow.
* added comments to demo scripts.
* fixed a bug in the quadratic operator  (a beta tester alerted us to this).
* means vectors are now always 1xn matrices (beta tester request).
* added an option "compute" to mxEval() to pre-compute matrix expressions without going to the optimizer.
* Matrix algebra conformability is now tested in R at the beginning of each mxRun().
* named entities (i.e. mxMatrices, mxAlgebras, etc.) can no longer have the same name as the label of a free parameter.  (This seems obscure, but you will like what we do with it in the next version!)
* can use options(mxByrow=TRUE) in the R global options if you always read your matrices in with the byrow=TRUE argument.  Saves some typing.  (beta tester request)
* fixed the standard error estimates summary.
* added mxVersion() function to return the version number (as a string).

Release 0.1.2-708 (Aug 14, 2009)
================================
* Added R help documentation for omxCheckCloseEnough(), omxCheckWithinPercentError(), 
	omxCheckTrue(), omxCheckEquals(), and omxCheckSetEquals()

* (mxMatrix) Fixed a bug in construction of symmetric matrixes.
    - now supports lower, standardized, and subdiagonal matrices.

Release 0.1 (Aug 03, 2009)
==========================
* (mxEvaluate) mxEvaluate translates MxMatrix references, MxAlgebra references,  MxObjectiveFunction references, and label references.

* (mxOptions) added 'reset' argument to mxOptions()

* (mxPath) renamed 'start' argument of mxPath() to 'values'
    - renamed 'name' argument of mxPath() to 'labels'
    - renamed 'boundMin' argument of mxPath() to 'lbound'
    - renamed 'boundMax' argument of mxPath() to 'ubound'
    - eliminated 'ciLower' argument of mxPath()
    - eliminated 'ciUpper' argument of mxPath()
    - eliminated 'description' argument of mxPath()

* (dimnames) implemented dimnames(x) for MxMatrix objects
    - implemented dimnames(x) <- value for MxMatrix objects
    - implemented dimnames(x) for MxAlgebra objects
    - implemented dimnames(x) <- value for MxAlgebra objects

* (mxMatrix) added 'dimnames' argument to mxMatrix()

* (mxData) renamed 'vector' argument of mxData() to 'means'
