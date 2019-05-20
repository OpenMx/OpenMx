# Where to find news

OpenMx developers, being lazy and incorrigible, often forget to update the NEWS file. To learn about new and exciting features, please visit https://openmx.ssri.psu.edu/

# OpenMx 2.13.0 FUTURE 2019 (R 3.6.0))
* NEW!: SEs for models with constraints! (used not to be calculated. Big win for twin models)
* IMPROVED: `omxSetParameters` warns if you don't ask to do anything.
* IMPROVED: `mxCheckIdentification` works with models containing constraints
* IMPROVED: `mxPath` more helpful messages for errors
* IMPROVED: `mxModel` error when path added to model without type = "RAM"
* IMPROVED: `summary` better at not counting redundant constraints.
* IMPROVED: progress reporting, inc. for mxTryHard.
* IMPROVED: github and CRAN README.md
* CHANGED: CXX14 (not CXX11)


# OpenMx 2.12.2 Feb 2019 (R 3.5.2))
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

# OpenMx 2.12.1 (Jan 20 2019 (R 3.5.2))
* Bug fix to definition variable handling in multilevel models (v2.11.5-2-g7ef03e3fe)

# OpenMx 2.11.4 (September 24 2018 (R 3.5.1))
* PARTYTIME:  Appears to be passing compiling for MacOS on CRAN !! 
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

# OpenMx 2.9.9 (R 3.5.0) Summer 2018
* No change for users

# OpenMx CRAN 2.9.6 (R 3.4.4)
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
	* If using 2.7.x (older), set expectation$.useSufficientSets = FALSE to avoid athis bug.

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
	* nb: multithreading supported on Linux and the OpenMx Team's build for Mac OS X.
* NEW: SLSQP multithreading to evaluate the numerical gradient.
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
