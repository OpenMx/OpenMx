.. _ordinal-specification-path:

Ordinal and Joint Ordinal-Continuous Model Specification
========================================================

This chapter deals with the specification of models that are either fit exclusively to ordinal variables or to a mix of ordinal and continuous variables. It extends the continuous data common factor model found in previous chapters to ordinal data.

The examples for this chapter can be found in the following file:

* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/OneFactorOrdinal_PathRaw.R
* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/OneFactorJoint_PathRaw.R

The continuous version of this model for raw data can be found the previous demos here:

* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/OneFactorModel_PathRaw.R

Ordinal Data Modeling
---------------------

OpenMx models ordinal data under a threshold model. A continuous normal distribution is assumed to underly every ordinal variable. These latent continuous distributions are only observed as being above or below a threshold, where there is one fewer threshold than observed categories in the data. For example, consider a variable with three ordered categories indicated by the values zero, one and two. Under this approach, this variable is assumed to follow a normal distribution that is partitioned or cut by two thresholds: individuals with underlying scores below the first threshold have an observed value of zero, individuals with latent scores between the thresholds are observed with values of one, and individuals with underlying scores give observed values of two.

.. image:: graph/ThresholdModel.png
	:height: 2in
	
Each threshold may be freely estimated or assigned as a fixed parameter, depending on the desired model. In addition to the thresholds, ordinal variables still have a mean and variance that describes the parameters of the underlying continuous distribution. However, this underlying distribution must be scaled by fixing at least two parameters to identify the model. One method of identification fixes the mean and variance to specific values, most commonly to a standard normal distribution with a mean of zero and a variance of one. A variation on this method fixes the residual variance of the categorical variable to one, which is often easier to specify. Alternatively, categorical variables may be identified by fixing two thresholds to non-equivalent constant values. These methods will differ in the scale assigned to the ordinal variables (and thus, the scale of the parameters estimated from them), but all identify the same model and should provide equally valid results.

OpenMx allows for the inclusion of continuous and ordinal variables in the same model, as well as models with only continuous or only ordinal variables. Any number of continuous variables may be included in an OpenMx model; however, maximum likelihood estimation for ordinal data must be limited to twenty ordinal variables regardless of the number of continuous variables. Further technical details on ordinal and joint continuous-ordinal optimization are contained at the end of this chapter.

Data Specification
^^^^^^^^^^^^^^^^^^

To use ordinal variables in OpenMx, users must identify ordinal variables by specifying those variables as ordered factors in the included data. Ordinal models can only be fit to raw data; if data is described as a covariance or other moment matrix, then the categorical nature of the data was already modeled to generate that moment matrix. Ordinal variables must be defined as specific columns in an R data frame.

Factors are a type of variable included in an R data frame. Unlike numeric or continuous variables, which must include only numeric and missing values, observed values for factors are treated as character strings. All factors contain a ``levels`` argument, which lists the possible values for a factor. Ordered factors contain information about the ordering of possible levels. Both R and OpenMx have tools for manipulating factors in data frames. The R functions ``factor()`` and ``as.factor()`` (and companions ``ordered()`` and ``as.ordered()``) can be used to specify ordered factors. OpenMx includes a helper function ``mxFactor()`` which more directly prepares ordinal variables as ordered factors in preparation for inclusion in OpenMx models. The code below demonstrates the ``mxFactor()`` function, replacing the variable *z1* that was initially read as a continuous variable and treating it as an ordinal variable with two levels. This process is repeated for *z2* (two levels) and *z3* (three levels).

.. cssclass:: input
..

.. code-block:: r

    data(myFADataRaw)

    oneFactorOrd <- myFADataRaw[,c("z1", "z2", "z3")]

    oneFactorOrd$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
    oneFactorOrd$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
    oneFactorOrd$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))

Threshold Specification
^^^^^^^^^^^^^^^^^^^^^^^

Just as covariances and means are included in models using ``mxPath`` when ``type='RAM'`` is enabled, thresholds may be included in models using the ``mxThreshold`` function. This function creates a list of thresholds to be added to your model, just as ``mxPath`` creates a list of paths. As an example, the data prep example above includes two binary variables (*z1* and *z2*) and one variable with three categories (*z3*). This means that models fit to this data should contain thresholds for three variables (for *z1*, *z2* and *z3*). This can be done with three separate calls to the ``mxThreshold`` function, as shown here.

.. cssclass:: input
..

.. code-block:: r

    mxThreshold(vars="z1", nThresh=1, free=TRUE, values=-1)
    mxThreshold(vars="z2", nThresh=1, free=TRUE, values=0)
    mxThreshold(vars="z3", nThresh=2, free=TRUE, values=c(-.5, 1.2))

The ``mxThreshold`` function first requires a variable to assign thresholds to, as well as a number of thresholds. In the first use of ``mxThreshold`` above, those are specified using the ``vars`` and ``nThresh`` arguments. The remaining arguments match those used by ``mxPath``: threshold parameters should be designated as ``free``, be given starting ``values``, and optionally given ``labels`` and boundaries (``lbound`` and ``ubound``). 
	
In this example, variables 'z1' and 'z2' are binary, with a single freely estimated threshold for each variable with starting values of -1 and 0, respectively. The meaning of these thresholds will depend on the mean and variance of these variables; as we are freely estimating thresholds for binary variables, the mean and variances of these variables should be constrained to fixed values. The third function call represents variable 'z3', which contains two thresholds and thus three categories. These two thresholds are assigned free parameters with staring values of -0.5 and 1.2, and the mean and variance of this variable should also be constrained to fixed values for identification. For variables with multiple thresholds, starting values should be monotonically increasing in each column such that the first column represents the first threshold and lowest value and the last column represents the last threshold and highest value.
	
Alternatively, ``mxThreshold`` can be used to specify thresholds for multiple variables at once. In the code below, ``mxThreshold`` is used to specify thresholds for all variables simultaneously. First, the ``vars`` argument contains a vector of variable names for which thresholds should be specified. The ``nThresh`` argument then specifies how many thresholds should be assigned to each variable: 1 each for *z1* and *z2*, and two for *z3*. The ``free`` argument states that all specified thresholds are to be freely estimated (the one value is repeated for all four thresholds). Finally, starting values are given using the ``values`` argument: -1 for *z1*, 0 for *z2*, and -.5 and 1.2 for *z3*.

.. cssclass:: input
..

.. code-block:: r

    mxThreshold(vars=c("z1","z2","z3"), nThresh=c(1,1,2), free=TRUE, values=c(-1,0,-.5,1.2) )

There are a few common errors regarding the use of thresholds in OpenMx. First, threshold values within each variable must be strictly increasing, such that the value in any element of the threshold matrix must be greater than all values above it in that column. In the above example, the second threshold for *z3* is set at 1.2, above the value of -.5 for the first threshold. OpenMx will return an error when your thresholds are not strictly increasing. There are no restrictions on values across variables: the second threshold for *z3* could be below all thresholds for *z1* and *z2* provided it exceeded the value for the first *z3* threshold. Second, the variables in your model that are assigned thresholds must match ordinal factors in the data. Additionally, free parameters should only be included for thresholds present in your data: including a second freely estimated threshold for *z1* or *z2* in this example would not directly impede model estimation, but would remain at its starting value and count as a free parameter for the purposes of calculating fit statistics.

It is also important to remember that specifying thresholds is not sufficient to get an ordinal data model to run. In addition, the scale of each ordinal variable must be identified just like the scale of a latent variable. The most common method for this involves constraining a ordinal item's mean to zero and either its total or residual variance to a constant value (i.e., one). For variables with two or more thresholds, ordinal variables may also be identified by constraining two thresholds to fixed values. Models that don't identify the scale of their ordinal variables should not converge.

Thresholds may also be expressed in matrix form. This is described in more detail in the matrix version of this chapter.

Users of original or ''classic'' Mx may recall specifying thresholds not in absolute terms, but as deviations. This method estimated the difference between each threshold for a variable and the previous one, which ensured that thresholds were in the correct order (i.e., that the second threshold for a variable was not lower than the first). While users may employ this method using ``mxAlgebra`` as it suits them, OpenMx does not require this technique. Simply specifying a thresholds matrix is typically sufficient to keep thresholds in proper order.

Including Thresholds in Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you use ``mxThreshold`` to specify thresholds, there is nothing left to do prior to running your model. However, if you manually create a threshold matrix, you must also specify the name of this matrix in your expectation function. This is described in more detail in the matrix version of this chapter.

Common Factor Model 
-------------------

All of the raw data examples through the documentation may be converted to ordinal examples by the inclusion of ordinal data, the specification of a threshold matrix and inclusion of that threshold matrix in the objective function. 

Ordinal Data
^^^^^^^^^^^^

The following example is a version of the continuous data common factor model referenced at the beginning of this chapter. Aside from replacing the continuous variables ``x1-x6`` with the ordinal variables ``z1-z3``, the code below simply incorporates the steps referenced above into the existing example. Data preparation occurs first, with the added ``mxFactor`` statements to identify ordinal variables and their ordered levels.

.. cssclass:: input
..

.. code-block:: r

    require(OpenMx)

    data(myFADataRaw)
    oneFactorOrd <- myFADataRaw[,c("z1","z2","z3")]

    oneFactorOrd$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0,1))
    oneFactorOrd$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0,1))
    oneFactorOrd$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0,1,2))

Model specification can be achieved by appending the above threshold matrix and objective function to either the path or matrix common factor examples. The path example below has been altered by changing the variable names from ``x1-x6`` to ``z1-z3``, adding the threshold matrix and objective function, and identifying the ordinal variables by constraining their means to be zero and their residual variances to be one.

.. cssclass:: input
..

.. code-block:: r

    dataRaw      <- mxData( observed=oneFactorOrd, type="raw" )
    # residual variances
    resVars      <- mxPath( from=c("z1","z2","z3"), arrows=2,
                            free=FALSE, values=c(1,1,1), labels=c("e1","e2","e3") )
    # latent variance
    latVar       <- mxPath( from="F1", arrows=2,
                            free=TRUE, values=1, labels ="varF1" )
    # factor loadings
    facLoads     <- mxPath( from="F1", to=c("z1","z2","z3"), arrows=1,
                            free=c(FALSE,TRUE,TRUE), values=1, labels=c("l1","l2","l3") )
    # means
    means        <- mxPath( from="one", to=c("z1","z2","z3","F1"), arrows=1,
                            free=FALSE, values=0, 
                            labels=c("meanz1","meanz2","meanz3","meanF") )
    # thresholds
    thresholds   <- mxThreshold( vars=c("z1","z2","z3"), nThresh=c(1,1,2), 
                            free=TRUE, values=c(-1,0,-.5,1.2) )
    oneFactorModel <- mxModel("Common Factor Model Path Specification", type="RAM",
                            manifestVars=c("z1","z2","z3"), latentVars="F1",
                            dataRaw, resVars, latVar, facLoads, means, thresholds)

This model may then be optimized using the ``mxRun`` command.

.. cssclass:: input
..

.. code-block:: r

    oneFactorResults <- mxRun(oneFactorModel)

Joint Ordinal-Continuous Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Models with both continuous and ordinal variables may be specified just like any other ordinal data model. Threshold matrices in these models should contain columns only for the ordinal variables, and should contain column names to designate which variables are to be treated as ordinal. In the example below, the one factor model above is estimated with three continuous variables (``x1-x3``) and three ordinal variables (``z1-z3``).

.. cssclass:: input
..

.. code-block:: r

    require(OpenMx)

    oneFactorJoint <- myFADataRaw[,c("x1","x2","x3","z1","z2","z3")]

    oneFactorJoint$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0,1))
    oneFactorJoint$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0,1))
    oneFactorJoint$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0,1,2))

    dataRaw      <- mxData( observed=oneFactorJoint, type="raw" )
    # residual variances
    resVars      <- mxPath( from=c("x1","x2","x3","z1","z2","z3"), arrows=2,
                            free=c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
                            values=1, labels=c("e1","e2","e3","e4","e5","e6") )
    # latent variance
    latVar       <- mxPath( from="F1", arrows=2,
                            free=FALSE, values=1, labels ="varF1" )
    # factor loadings
    facLoads     <- mxPath( from="F1", to=c("x1","x2","x3","z1","z2","z3"), arrows=1,
                            free=TRUE, values=1, labels=c("l1","l2","l3","l4","l5","l6") )
    # means
    means        <- mxPath( from="one", to=c("x1","x2","x3","z1","z2","z3","F1"), arrows=1,
                            free=c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE), values=0,
                            labels=c("meanx1","meanx2","meanx3",
                                     "meanz1","meanz2","meanz3","meanF") )
    # thresholds
    thresholds   <- mxThreshold(vars=c("z1","z2","z3"), nThresh=c(1,1,2),
                            free=TRUE, values=c(-1,0,-.5,1.2) )
    oneFactorJointModel <- mxModel("Common Factor Model Path Specification", type="RAM",
                            manifestVars=c("x1","x2","x3","z1","z2","z3"), latentVars="F1",
                            dataRaw, resVars, latVar, facLoads, means, thresholds)

This model may then be optimized using the ``mxRun`` command.

.. cssclass:: input
..

.. code-block:: r

    oneFactorJointResults <- mxRun(oneFactorJointModel)

Technical Details
-----------------

Maximum likelihood estimation for ordinal variables is done by generating expected covariance and mean matrices for the latent continuous variables underlying the set of ordinal variables, then integrating the multivariate normal distribution defined by those covariances and means. The likelihood for each row of the data is defined as the multivariate integral of the expected distribution over the interval defined by the thresholds bordering that row's data. OpenMx uses Alan Genz's SADMVN routine for multivariate normal integration (see http://www.math.wsu.edu/faculty/genz/software/software.html for more information). 

When continuous variables are present, OpenMx utilizes a block decomposition to separate the continuous and ordinal covariance matrices for FIML. The likelihood of the continuous variables is calculated normally.  The effects of the point estimates of the continuous variables is projected out of the expected covariance matrix of the ordinal data. The likelihood of the ordinal data is defined as the multivariate integral over the distribution defined by the resulting ordinal covariance matrix.
