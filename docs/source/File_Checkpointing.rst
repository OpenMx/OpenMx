.. _file-checkpointing:

File Checkpointing
==================

This section will cover how to periodically save the state of the optimizer to a file on disk.  In the event of a system crash, the checkpoint file can be loaded back into the model when the system has been restored.  The checkpoint file can also be used as a log to trace the path of optimization. The simplest form of file checkpointing is to use the argument ``checkpoint = TRUE`` in the call to the ``mxRun()`` function.

.. code-block:: r

	require(OpenMx)

	data(demoOneFactor)
	manifestVars <- names(demoOneFactor)

	factorModel <- mxModel("One Factor",
	    mxMatrix(type="Full", nrow=5, ncol=1, values=0.2, free=TRUE, name="A", labels=letters[1:5]),
	    mxMatrix(type="Symm", nrow=1, ncol=1, values=1, free=FALSE, name="L"),
	    mxMatrix(type="Diag", nrow=5, ncol=5, values=1, free=TRUE, name="U"),
	    mxAlgebra(expression=A %*% L %*% t(A) + U, name="R"),
	    mxMLObjective(covariance="R", dimnames=manifestVars),
	    mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
	)

	factorFit <- mxRun(factorModel, checkpoint = TRUE)

With no extra options, a checkpoint file will be created in the current working directory with the filename: "<modelname>.omx". The checkpoint file is a data.frame object such that each row contains all the values of the free parameters at a particular instance in time. By default, a row is added to the file every 10 minutes.  The ``mxOption()`` function can be used to set the directory of the checkpoint file, an optional prefix to the checkpoint filename, the decision to save based on minutes or optimizer iterations, and the checkpoint interval in either minutes or optimizer iterations. Below is an example that modifies some of the checkpoint options:

.. code-block:: r

	directory <- tempdir()

	factorModel <- mxOption(factorModel, "Checkpoint Directory", directory)
	factorModel <- mxOption(factorModel, "Checkpoint Units", "iterations")
	factorModel <- mxOption(factorModel, "Checkpoint Count", 10)

After a checkpoint file has been created, it can be loaded into a MxModel object using the ``mxRestore()`` function.  It is necessary to specify the checkpoint directory and checkpoint filename prefix if they were declared when the checkpoint file was created:

.. code-block:: r
	
	factorFit <- mxRun(factorModel, checkpoint = TRUE)	
	factorRestore <- mxRestore(factorModel, chkpt.directory = directory)

The checkpoint directory will extend to independent submodels in a collection of models.  Each independent submodel will be saved in a separate file.  See ``?mxOption`` or ``getOption('mxOptions')`` for a list of options that modify the behavior of file checkpointing.  See ``?mxRestore`` for more information on restoring a checkpoint file.
