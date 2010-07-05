.. _multicore-execution:

Multicore Execution
===================

This section will cover how take advantage of multiple cores on your machine.  To use the multicore mode in OpenMx, you must declare independent submodels. A model is declared independent by using the argument 'independent=TRUE' in the ``mxModel()`` function. An independent model and all of its dependent children are executed in a separate optimization environment. An independent model shares **no** free parameters with either its sibling models or its parent model. An independent model may **not** refer to matrices or algebras in either its sibling models or its parent model. A parent model may access the final results of optimization from an independent child model. 

To use the snowfall library, you must start your R environment with the following commands:

.. code-block:: r

	library(OpenMx)
	library(snowfall)
	sfInit(parallel = TRUE, cpus = ???)
	sfLibrary(OpenMx)

``sfInit`` will initialize the snowfall cluster. You must specify either the number of CPUs on your machine or the cluster environment (see snowfall package documentation). ``sfLibrary`` exports the OpenMx library to the client nodes in the cluster. At the end of your script, use the command:

.. code-block:: r

	sfStop()

To Improve Performance
----------------------

Any sequential portions of your script will quickly become the performance bottleneck (Amdahl's Law). Avoid iteration over large data structures. Use the functions ``omxLapply()`` and ``omxSapply()`` instead of iteration. These two functions invoke the snowfall ``sfLapply()`` and ``sfSapply()`` functions if the snowfall library has been loaded. Otherwise they invoke the sequential functions ``lapply()`` and ``sapply()``. To hunt for bottelnecks in your script, run your script with multicore settings enabled and use Rprof to profile a reasonable size test case. Ignore the calls to sfLapply() and sfSapply() in the results of profiling. Any other time-consuming calls represent potential sequential bottlenecks.

Some of the functions provided by the OpenMx library can be bottlenecks. Iterative use of the ``mxModel()`` function in order to add submodels can be time consuming. Use the following unsafe idiom to improve performance:

.. code-block:: r

	# generate a list of independent submodels
	submodels <- # omxLapply(...)
	names(submodels) <- omxExtractNames(submodels)
	topModel@submodels <- submodels

An Example
----------

The following script can be found with ``demo(BootstrapParallel)``

.. code-block:: r

	# parameters for the simulation: lambda = factor loadings,
	# specifics = specific variances
	lambda <- matrix(c(.8, .5, .7, 0), 4, 1)
	nObs <- 500
	nReps <- 10
	nVar <- nrow(lambda)
	specifics <- diag(nVar)
	chl <- chol(lambda %*% t(lambda) + specifics)

	# indices for parameters and hessian estimate in results
	pStrt <- 3
	pEnd <- pStrt + 2*nVar - 1
	hStrt <- pEnd + 1
	hEnd <- hStrt + 2*nVar - 1

	# dimension names for OpenMx
	dn <- list()
	dn[[1]] <- paste("Var", 1:4, sep="")
	dn[[2]] <- dn[[1]]

	# function to get a covariance matrix
	randomCov <- function(nObs, nVar, chl, dn) {
		x <- matrix(rnorm(nObs*nVar), nObs, nVar)
		x <- x %*% chl
		thisCov <- cov(x)
		dimnames(thisCov) <- dn
		return(thisCov)  
	}

	createNewModel <- function(index, prefix, model) {
		modelname <- paste(prefix, index, sep='')
		model@data@observed <- randomCov(nObs, nVar, chl, dn)
		model@name <- modelname
		return(model)
	}

	getStats <- function(model) {
		retval <- c(model@output$status[[1]],
			max(abs(model@output$gradient)),
			model@output$estimate,
			sqrt(diag(solve(model@output$hessian))))
		return(retval)
	}


	# initialize obsCov for MxModel
	obsCov <- randomCov(nObs, nVar, chl, dn)

	# results matrix: get results for each simulation
	results <- matrix(0, nReps, hEnd)
	dnr <- c("inform", "maxAbsG", paste("lambda", 1:nVar, sep=""),
		paste("specifics", 1:nVar, sep=""),
		paste("hessLambda", 1:nVar, sep=""),
		paste("hessSpecifics", 1:nVar, sep=""))
	dimnames(results)[[2]] <- dnr

	# instantiate MxModel
	template <- mxModel(name="stErrSim",
                       mxMatrix(name="lambda", type="Full", nrow=4, ncol=1,
                                free=TRUE, values=c(.8, .5, .7, 0)),
                       mxMatrix(name="specifics", type="Diag", nrow=4,
                                free=TRUE, values=rep(1, 4)),
                       mxAlgebra(lambda %*% t(lambda) + specifics,
                                 name="preCov", dimnames=dn),
                       mxData(observed=obsCov, type="cov", numObs=nObs),
                       mxMLObjective(covariance='preCov'),
                       independent = TRUE)

	topModel <- mxModel(name = 'container')

	submodels <- lapply(1:nReps, createNewModel, 'stErrSim', template)

	names(submodels) <- omxExtractNames(submodels)
	topModel@submodels <- submodels

	modelResults <- mxRun(topModel, silent=TRUE, suppressWarnings=TRUE)

	results <- t(omxSapply(modelResults@submodels, getStats))

	# get rid of bad covergence results
	results2 <- data.frame(results[which(results[,1] <= 1),])

	# summarize the results
	means <- mean(results2)
	stdevs <- sd(results2)
	sumResults <- data.frame(matrix(dnr[pStrt:pEnd], 2*nVar, 1,
                                dimnames=list(NULL, "Parameter")))
	sumResults$mean <- means[pStrt:pEnd]
	sumResults$obsStDev <- stdevs[pStrt:pEnd]
	sumResults$meanHessEst <- means[hStrt:hEnd]
	sumResults$sqrt2meanHessEst <- sqrt(2) * sumResults$meanHessEst

	# print results
	print(sumResults)

