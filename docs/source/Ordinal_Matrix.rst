.. _ordinal-specification-matrix:

Ordinal and Joint Ordinal-Continuous Model Specification
========================================================

This chapter deals with the specification of models that are either fit exclusively to ordinal variables or to a mix of ordinal and continuous variables. It extends the continuous data common factor model found in previous chapters to ordinal data.

The examples for this chapter can be found in the following files:

* http://openmx.psyc.virginia.edu/svn/trunk/demo/OneFactorOrdinal_MatrixRaw.R
* http://openmx.psyc.virginia.edu/svn/trunk/demo/OneFactorOrdinal01_MatrixRaw.R

The continuous versions of these models for raw data can be found the previous demos here:

* http://openmx.psyc.virginia.edu/svn/trunk/demo/OneFactorModel_MatrixRaw.R
* http://openmx.psyc.virginia.edu/svn/trunk/demo/OneFactorModel_PathRaw.R

Ordinal Data
------------

OpenMx models ordinal data under a threshold model. A continuous normal distribution is assumed to underly every ordinal variable. These latent continuous distributions are only observed as being above or below a threshold, where there is one fewer threshold than observed categories in the data. For example, consider a variable with three ordered categories indicated by the values zero, one and two. Under this approach, this variable is assumed to follow a normal distribution that is partitioned or cut by two thresholds: individuals with underlying scores below the first threshold have an observed value of zero, individuals with latent scores between the thresholds are observed with values of one, and individuals with underlying scores give observed values of two.

.. image:: graph/thresh.png
    :height: 2in

Each threshold may be freely estimated or assigned as a fixed parameter, depending on the desired model. In addition to the thresholds, ordinal variables still have a mean and variance that describes the parameters of the underlying continuous distribution. However, this underlying distribution must be scaled by fixing at least two parameters to identify the model. One method of identification fixes the mean and variance to specific values, most commonly to a standard normal distribution with a mean of zero and a variance of one. A variation on this method fixes the residual variance of the categorical variable to one, which is often easier to specify. Alternatively, categorical variables may be identified by fixing two thresholds to non-equivalent constant values. These methods will differ in the scale assigned to the ordinal variables (and thus, the scale of the parameters estimated from them), but all identify the same model and should provide equally valid results.

OpenMx allows for the inclusion of continuous and ordinal variables in the same model, as well as models with only continuous or only ordinal variables. Any number of continuous variables may be included in an OpenMx model; however, maximum likelihood estimation for ordinal data must be limited to twenty ordinal variables regardless of the number of continuous variables. Further technical details on ordinal and joint continuous-ordinal optimization are contained at the end of this chapter.

Specifying Data for Ordinal Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use ordinal variables in OpenMx, users must identify ordinal variables by specifying those variables as ordered factors in the included data. Ordinal models can only be fit to raw data; if data is described as a covariance or other moment matrix, then the categorical nature of the data was already modeled to generate that moment matrix. Ordinal variables must be defined as specific columns in an R data frame.

Factors are a type of variable included in an R data frame. Unlike numeric or continuous variables, which must include only numeric and missing values, observed values for factors are treated as character strings. All factors contain a ``levels`` argument, which lists the possible values for a factor. Ordered factors contain information about the ordering of possible levels. Both R and OpenMx have tools for manipulating factors in data frames. The R functions ``factor()`` and ``as.factor()`` (and companions ``ordered()`` and ``as.ordered()``) can be used to specify ordered factors. OpenMx includes a helper function ``mxFactor()`` which more directly prepares ordinal variables as ordered factors in preparation for inclusion in OpenMx models. The code below demonstrates the ``mxFactor()`` function, replacing the variable *z1* that was initially read as a continuous variable and treating it as an ordinal variable with two levels. This process is repeated for *z2* (two levels) and *z3* (three levels).

.. code-block:: r

    data(myFADataRaw)

    oneFactorOrd <- myFADataRaw[,c("z1", "z2", "z3")]

    oneFactorOrd$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
    oneFactorOrd$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
    oneFactorOrd$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))

Specifying Threshold Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Just as covariances and means are included in models by specifying matrices and algebras, thresholds may be included in models as threshold matrices. These matrices can be of user-specified type, though most will be of type ``Full``. The columns of this matrix should correspond to the ordinal variables in your dataset, with the column names of this matrix corresponding to variables in your data. This assignment can be done either with the ``dimnames`` argument to ``mxMatrix``, or by using the ``threshnames`` argument in your expectation function the same way ``dimnames`` arguments are used. The rows of your threshold matrix should correspond to the ordered thresholds for each variable, such that the first row is the lowest threshold for each variable, the second row is the next threshold (provided one or more of your variables have two thresholds), and so on for the maximum number of thresholds you have in your data. Rows of the threshold matrix beyond the number of thresholds in a particular variable should be fixed parameters with starting values of ``NA``.

As an example, the data prep example above includes two binary variables (*z1* and *z2*) and one variable with three categories (*z3*). This means that the threshold matrix for models fit to this data should contain three columns (for *z1*, *z2* and *z3*) and two rows, as the variable *z3* requires two thresholds. The code below specifies a 2 x 3 ``Full`` matrix with free parameters for one threshold for *z1*, one threshold for *z2* and two thresholds for *z3*.

.. code-block:: r

    thresh       <- mxMatrix( type="Full", nrow=2, ncol=3,
                              free=c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE), 
                              values=c(-1,0,-.5,NA,NA,1.2), byrow=TRUE, name="thresh" )

There are a few common errors regarding the use of thresholds in OpenMx. First, threshold values within each row must be strictly increasing, such that the value in any element of the threshold matrix must be greater than all values above it in that column. In the above example, the second threshold for *z3* is set at 1.2, above the value of -.5 for the first threshold. OpenMx will return an error when your thresholds are not strictly increasing. There are no restrictions on values across columns or variables: the second threshold for *z3* could be below all thresholds for *z1* and *z2* provided it exceeded the value for the first *z3* threshold. Second, the dimnames of the threshold matrix must match ordinal factors in the data. Additionally, free parameters should only be included for thresholds present in your data: including a second freely estimated threshold for *z1* or *z2* in this example would not directly impede model estimation, but would remain at its starting value and count as a free parameter for the purposes of calculating fit statistics.

It is also important to remember that specifying a threshold matrix is not sufficient to get an ordinal data model to run. In addition, the scale of each ordinal variable must be identified just like the scale of a latent variable. The most common method for this involves constraining a ordinal item's mean to zero and either its total or residual variance to a constant value (i.e., one). For variables with two or more thresholds, ordinal variables may also be identified by constraining two thresholds to fixed values. Models that don't identify the scale of their ordinal variables should not converge.

While thresholds can't be expressed as paths between variables like other parts of the model, OpenMx supports a path-like interface called ``mxThreshold`` as of version 2.0. This function is described in more detail in the ordinal data version of this chapter and the ``mxThreshold`` help file.

Users of original or ''classic'' Mx may recall specifying thresholds not in absolute terms, but as deviations. This method estimated the difference between each threshold for a variable and the previous one, which ensured that thresholds were in the correct order (i.e., that the second threshold for a variable was not lower than the first). While users may employ this method using ``mxAlgebra`` as it suits them, OpenMx does not require this technique. Simply specifying a thresholds matrix is typically sufficient to keep thresholds in proper order.
	
Including Thresholds in Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, the threshold matrix must be identified as such in the expectation function in the same way that other matrices are identified as means or covariance matrices. Both the ``mxExpectationNormal`` and ``mxExpectationRAM`` contain a ``thresholds`` argument, which takes the name of the matrix or algebra to be used as the threshold matrix for a given analysis. Although specifying ``type='RAM'`` generates a RAM expectation function, this expectation function must be replaced by one with a specified thresholds matrix.

You must specify ``dimnames`` (dimension names) for your thresholds matrix that correspond to the ordered factors in the data you wish to analyze. This may be done in either of two ways, both of which correspond to specifying dimnames for other OpenMx matrices. One method is to use the ``threshnames`` argument in the ``mxExpectationNormal`` or ``mxExpectationRAM`` functions, which specifies which variables are in a threshold matrix in the same way the ``dimnames`` argument specifies which variables are in the rest of the model. Another method is to specify dimnames for each matrix using the ``dimnames`` argument in the ``mxMatrix`` function. Either method may be used, but it is important to use the same method for all matrices in a given model (either using expectation function arguments ``dimnames`` and ``threshnames`` or supplying ``dimnames`` for all ``mxMatrix`` objects manually). Expectation function arguments ``dimnames`` and ``threshnames`` supersede the matrix ``dimname`` arguments, and ``threshnames`` will take the value of the ``dimnames`` if both ``dimnames`` and ``thresholds`` are specified but ``threshnames`` is omitted. 

The code below specifies an ``mxExpectationRAM`` to include a thresholds matrix named ``"thresh"``. When models are built using ``type='RAM'``, the ``dimnames`` argument may be omitted, as the requisite dimnames for the ``A``, ``S``, ``F`` and ``M`` matrices are generated from the ``manifestVars`` and ``latentVars`` lists. However, the dimnames for the threshold matrix should be included using the ``dimnames`` argument in ``mxMatrix``.

.. code-block:: r

	mxExpectationRAM(A="A", S="S", F="F", M="M", thresholds="thresh")

Common Factor Model for Ordinal Data
-_----------------------------------

All of the raw data examples through the documentation may be converted to ordinal examples by the inclusion of ordinal data, the specification of a threshold matrix and inclusion of that threshold matrix in the expectation function. The following example is a version of the continuous data common factor model referenced at the beginning of this chapter. Aside from replacing the continuous variables ``x1-x6`` with the ordinal variables ``z1-z3``, the code below simply incorporates the steps referenced above into the existing example. Data preparation occurs first, with the added ``mxFactor`` statements to identify ordinal variables and their ordered levels.

.. code-block:: r

    require(OpenMx)

    data(myFADataRaw)

    oneFactorOrd <- myFADataRaw[,c("z1", "z2", "z3")]

    oneFactorOrd$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
    oneFactorOrd$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
    oneFactorOrd$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))

Model specification can be achieved by appending the above threshold matrix and expectation function to either the path or matrix common factor examples. The path example below has been altered by changing the variable names from ``x1-x6`` to ``z1-z3``, adding the threshold matrix and expectation function, and identifying the ordinal variables by constraining their means to be zero and their residual variances to be one.

.. code-block:: r

    dataRaw      <- mxData(oneFactorOrd, type="raw"),
    # asymmetric paths
    matrA        <- mxMatrix( type="Full", nrow=4, ncol=4,
                              free=c(F,F,F,T,
                                     F,F,F,T,
                                     F,F,F,T,
                                     F,F,F,F),
                              values=c(0,0,0,1,
                                       0,0,0,1,
                                       0,0,0,1,
                                       0,0,0,0),
                              labels=c(NA,NA,NA,"l1",
                                       NA,NA,NA,"l2",
                                       NA,NA,NA,"l3",
                                       NA,NA,NA,NA),
                              byrow=TRUE, name="A" )
    # symmetric paths
    matrS        <- mxMatrix( type="Symm", nrow=4, ncol=4, 
                              free=FALSE, 
                              values=diag(4)
                              labels=c("e1", NA, NA,  NA,
                                        NA,"e2", NA,  NA,
                                        NA,  NA,"e3", NA,
                                        NA,  NA, NA, "varF1"),
                              byrow=TRUE, name="S" )
    # filter matrix
    matrF        <- mxMatrix( type="Full", nrow=3, ncol=4,
                              free=FALSE, values=c(1,0,0,0,  0,1,0,0,  0,0,1,0),
                              byrow=TRUE, name="F" )
    # means
    matrM        <- mxMatrix( type="Full", nrow=1, ncol=4,
                              free=FALSE, values=0,
                              labels=c("meanz1","meanz2","meanz3",NA), name="M" )
    thresh       <- mxMatrix( type="Full", nrow=2, ncol=3,
                              free=c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE), 
                              values=c(-1,0,-.5,NA,NA,1.2), byrow=TRUE, name="thresh" )
    exp          <- mxExpectationRAM("A","S","F","M", dimnames=c("z1","z2","z3","F1"), 
                              thresholds="thresh", threshnames=c("z1","z2","z3")),
    funML        <- mxFitFunctionML()

    oneFactorModel <- mxModel("Common Factor Model Matrix Specification", 
                               dataRaw, matrA, matrS, matrF, matrM, thresh, exp, funML)
                           
This model may then be optimized using the ``mxRun`` command.

.. code-block:: r

    oneFactorResults <- mxRun(oneFactorModel)

Common Factor Model for Joint Ordinal-Continuous Data
-----------------------------------------------------

Models with both continuous and ordinal variables may be specified just like any other ordinal data model. Threshold matrices in these models should contain columns only for the ordinal variables, and should contain column names to designate which variables are to be treated as ordinal. In the example below, the one factor model above is estimated with three continuous variables (``x1-x3``) and three ordinal variables (``z1-z3``).

.. code-block:: r

    require(OpenMx)

    oneFactorJoint <- myFADataRaw[,c("x1", "x2", "x3", "z1", "z2", "z3")]

    oneFactorJoint$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
    oneFactorJoint$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
    oneFactorJoint$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))

    dataRaw      <- mxData(observed=oneFactorJoint, type="raw")
    # asymmetric paths
    matrA        <- mxMatrix( type="Full", nrow=7, ncol=7,
                              free=c(rep(c(F,F,F,F,F,F,T),6),rep(F,7)),
                              values=c(rep(c(0,0,0,0,0,0,1),6),rep(F,7)),
                              labels=rbind(cbind(matrix(NA,6,6),matrix(paste("l",1:6,sep=""),6,1)),
                               matrix(NA,1,7)),
                              byrow=TRUE, name="A" )
    # symmetric paths
    labelsS      <- matrix(NA,7,7); diag(labelsS) <- c(paste("e",1:6,sep=""),"varF1")
    matrS        <- mxMatrix( type="Symm", nrow=7, ncol=7, 
                              free= rbind(cbind(matrix(as.logical(diag(3)),3,3),matrix(F,3,4)), 
                               matrix(F,4,7)),
                              values=diag(7), labels=labelsS, byrow=TRUE, name="S" )
    # filter matrix
    matrF        <- mxMatrix( type="Full", nrow=6, ncol=7,
                              free=FALSE, values=cbind(diag(6),matrix(0,6,1)),
                              byrow=TRUE, name="F" )
    # means
    matrM        <- mxMatrix( type="Full", nrow=1, ncol=7,
                              free=c(T,T,T,F,F,F,F), values=c(1,1,1,0,0,0,0),
                              labels=c("meanx1","meanx2","meanx3","meanz1","meanz2","meanz3",NA),
                              name="M" )
    thresh       <- mxMatrix( type="Full", nrow=2, ncol=3,
                              free=c(TRUE,TRUE,TRUE,FALSE,FALSE,TRUE), 
                              values=c(-1,0,-.5,NA,NA,1.2), byrow=TRUE, name="thresh" )
    exp          <- mxExpectationRAM("A","S","F","M", 
                                     dimnames=c("x1","x2","x3","z1","z2","z3","F1"), 
                                     thresholds="thresh", threshnames=c("z1","z2","z3"))
    funML        <- mxFitFunctionML()

    oneFactorJointModel <- mxModel("Common Factor Model Matrix Specification", 
                                   dataRaw, matrA, matrS, matrF, matrM, thresh, exp, funML)
This model may then be optimized using the ``mxRun`` command.

.. code-block:: r

    oneFactorJointResults <- mxRun(oneFactorJointModel)

Technical Details
-----------------

Maximum likelihood estimation for ordinal variables by generating expected covariance and mean matrices for the latent continuous variables underlying the set of ordinal variables, then integrating the multivariate normal distribution defined by those covariances and means. The likelihood for each row of the data is defined as the multivariate integral of the expected distribution over the interval defined by the thresholds bordering that row's data. OpenMx uses Alan Genz's SADMVN routine for multivariate normal integration (see http://www.math.wsu.edu/faculty/genz/software/software.html for more information). 

When continuous variables are present, OpenMx utilizes a block decomposition to separate the continuous and ordinal covariance matrices for FIML. The likelihood of the continuous variables is calculated normally.  The effects of the point estimates of the continuous variables is projected out of the expected covariance matrix of the ordinal data. The likelihood of the ordinal data is defined as the multivariate integral over the distribution defined by the resulting ordinal covariance matrix.
