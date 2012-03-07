.. _multiplegroups-matrix-specification:

Multiple Groups, Matrix Specification
=====================================

An important aspect of structural equation modeling is the use of multiple groups to compare means and covariances structures between any two (or more) data groups, for example males and females, different ethnic groups, ages etc.  Other examples include groups which have different expected covariances matrices as a function of parameters in the model, and need to be evaluated together for the parameters to be identified.

The example includes the heterogeneity model as well as its submodel, the homogeneity model and is available in the following file:

* http://openmx.psyc.virginia.edu/svn/trunk/demo/BivariateHeterogeneity_MatrixRaw.R

A parallel version of this example, using paths specification of models rather than matrices, can be found here:

* http://openmx.psyc.virginia.edu/svn/trunk/demo/BivariateHeterogeneity_PathRaw.R


Heterogeneity Model
-------------------

We will start with a basic example here, building on modeling means and variances in a saturated model.  Assume we have two groups and we want to test whether they have the same mean and covariance structure.  

The path diagram of the heterogeneity model for a set of variables :math:`x` and :math:`y` are given below.

.. math::
..   :nowrap:
   
..   \begin{eqnarray*} 
..   x = \mu_{x1} + \sigma_{x1}
..   \end{eqnarray*}

.. image:: graph/HeterogeneityModel.png
    :height: 2.5in

Data
^^^^

For this example we simulated two datasets (``xy1`` and ``xy2``) each with zero means and unit variances, one with a correlation of 0.5, and the other with a correlation of 0.4 with 1000 subjects each.  We use the ``mvrnorm`` function in the ``MASS`` package, which takes three arguments: ``Sample Size``, ``Means``, ``Covariance Matrix``).  We check the means and covariance matrix in R and provide ``dimnames`` for the dataframe.  See attached R code for simulation and data summary.

.. code-block:: r

    #Simulate Data
    require(MASS)
    #group 1
    set.seed(200)
    rs=0.5
    xy1 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
    set.seed(200)
    #group 2
    rs=0.4
    xy2 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))

    #Print Descriptive Statistics
    selVars <- c('X','Y')
    summary(xy1)
    cov(xy1)
    dimnames(xy1) <- list(NULL, selVars)
    summary(xy2)
    cov(xy2)
    dimnames(xy2) <- list(NULL, selVars)
    
    
Model Specification
^^^^^^^^^^^^^^^^^^^

As before, we include the OpenMx package using a ``require`` statement.
We first fit a heterogeneity model, allowing differences in both the mean and covariance structure of the two groups.  As we are interested whether the two structures can be equated, we have to specify the models for the two groups, named ``group1`` and ``group2`` within another model, named ``bivHet``.  The structure of the job thus look as follows, with two ``mxModel`` commands as arguments of another ``mxModel`` command.  Note that ``mxModel`` commands are unlimited in the number of arguments.

.. code-block:: r

    require(OpenMx)

    bivHetModel <- mxModel("bivHet",
         mxModel("group1"),
         mxModel("group2"),
         mxAlgebra(group1.objective + group2.objective, name="minus2loglikelihood"),
         mxAlgebraObjective("minus2loglikelihood")
    )
     
For each of the groups, we fit a saturated model, using a Cholesky decomposition to generate the expected covariance matrix and a row vector for the expected means.  Note that we have specified different labels for all the free elements, in the two ``mxModel`` statements.  For more details, see example 1.

.. code-block:: r

    #Fit Heterogeneity Model
    bivHetModel <- mxModel("bivHet",
        mxModel("group1",
            mxMatrix(
                type="Lower", 
                nrow=2, 
                ncol=2, 
                free=T, 
                values=.5,
                labels=c("Ch11", "Ch21", "Ch31"),
                name="Chol1"
            ), 
            mxAlgebra(
                Chol1 %*% t(Chol1), 
                name="EC1" 
            ), 
            mxMatrix(
                type="Full", 
                nrow=1, 
                ncol=2, 
                free=T, 
                values=c(0,0), 
                labels=c("mX1", "mY1"), 
                name="EM1"
            ), 
            mxData(
                xy1, 
                type="raw"
            ), 
            mxFIMLObjective(
                covariance="EC1", 
                means="EM1",
                dimnames=selVars
            )
        ),
        mxModel("group2",
            mxMatrix(
                type="Lower", 
                nrow=2, 
                ncol=2, 
                free=T, 
                values=.5,
                labels=c("Ch12", "Ch22", "Ch32"),
                name="Chol2"
            ), 
            mxAlgebra(
                Chol2 %*% t(Chol2), 
                name="EC2"
            ), 
            mxMatrix(
                type="Full", 
                nrow=1, 
                ncol=2, 
                free=T, 
                values=c(0,0), 
                labels=c("mX2", "mY2"), 
                name="EM2"
            ), 
            mxData(
                xy2, 
                type="raw"
            ), 
            mxFIMLObjective(
                covariance="EC2", 
                means="EM2",
                dimnames=selVars
            )
        ),


We estimate five parameters (two means, two variances, one covariance) per group for a total of 10 free parameters.  We cut the ``Labels matrix:`` parts from the output generated with ``bivHetModel$group1@matrices`` and ``bivHetModel$group2@matrices``::

    in group1
        $S
                X      Y     
        X  "Ch11"     NA
        Y  "Ch21"  "Ch22" 

        $M
                X      Y    
        [1,] "mX1" "mY1"

    in group2
        $S
                X      Y     
        X  "Ch12"     NA
        Y  "Ch22" "Ch32" 

        $M
                X      Y    
        [1,] "mX2" "mY2"

To evaluate both models together, we use an ``mxAlgebra`` command that adds up the values of the objective functions of the two groups.  The objective function to be used here is the ``mxAlgebraObjective`` which uses as its argument the sum of the function values of the two groups, referred to by the name of the previously defined ``mxAlgebra`` object ``h12``.

.. code-block:: r

        mxAlgebra(
            group1.objective + group2.objective, 
            name="minus2loglikelihood"
        ),
        mxAlgebraObjective("minus2loglikelihood")
    )

Model Fitting
^^^^^^^^^^^^^

The ``mxRun`` command is required to actually evaluate the model.  Note that we have adopted the following notation of the objects.  The result of the ``mxModel`` command ends in "Model"; the result of the ``mxRun`` command ends in "Fit".  Of course, these are just suggested naming conventions.

.. code-block:: r

    bivHetFit <- mxRun(bivHetModel)

A variety of output can be printed.  We chose here to print the expected means and covariance matrices for the two groups and the likelihood of data given the model.  The ``mxEval`` command takes any R expression, followed by the fitted model name.  Given that the model ``bivHetFit`` included two models (group1 and group2), we need to use the two level names, i.e. ``group1.EM1`` to refer to the objects in the correct model.

.. code-block:: r

    EM1Het <- mxEval(group1.EM1, bivHetFit)
    EM2Het <- mxEval(group2.EM2, bivHetFit)
    EC1Het <- mxEval(group1.EC1, bivHetFit)
    EC2Het <- mxEval(group2.EC2, bivHetFit)
    LLHet <- mxEval(objective, bivHetFit)


Homogeneity Model: a Submodel
-----------------------------

Next, we fit a model in which the mean and covariance structure of the two groups are equated to one another, to test whether there are significant differences between the groups.  Rather than having to specify the entire model again, we copy the previous model ``bivHetModel`` into a new model ``bivHomModel`` to represent homogeneous structures.

.. code-block:: r

    #Fit Homogeneity Model
    bivHomModel <- bivHetModel

As elements in matrices can be equated by assigning the same label, we now have to equate the labels of the free parameters in group 1 to the labels of the corresponding elements in group 2.  This can be done by referring to the relevant matrices using the ``ModelName$MatrixName`` syntax, followed by ``@labels``.  Note that in the same way, one can refer to other arguments of the objects in the model.  Here we assign the labels from group1 to the labels of group2, separately for the Cholesky matrices used for the expected covariance matrices and for the expected means vectors.

.. code-block:: r

    bivHomModel$group2.Chol2@labels <- bivHomModel$group1.Chol1@labels
    bivHomModel$group2.EM2@labels <- bivHomModel$group1.EM1@labels

The specification for the submodel is reflected in the names of the labels which are now equal for the corresponding elements of the mean and covariance matrices, as below::

    in group1
        $S
                X      Y     
        X  "Ch11"     NA
        Y  "Ch21" "CH31" 

        $M
                X      Y    
        [1,] "mX1" "mY1"
    
    in group2
        $S
                X      Y     
        X  "Ch11"     NA
        Y  "Ch21" "Ch31" 

        $M
                X      Y     
        [1,] "mX1" "mY1"

We can produce similar output for the submodel, i.e. expected means and covariances and likelihood, the only difference in the code being the model name.  Note that as a result of equating the labels, the expected means and covariances of the two groups should be the same.

.. code-block:: r

    bivHomFit <- mxRun(bivHomModel)
        EM1Hom <- mxEval(group1.EM1, bivHomFit)
        EM2Hom <- mxEval(group2.EM2, bivHomFit)
        EC1Hom <- mxEval(group1.EC1, bivHomFit)
        EC2Hom <- mxEval(group2.EC2, bivHomFit)
        LLHom <- mxEval(objective, bivHomFit)

Finally, to evaluate which model fits the data best, we generate a likelihood ratio test as the difference between -2 times the log-likelihood of the homogeneity model and -2 times the log-likelihood of the heterogeneity model.  This statistic is asymptotically distributed as a Chi-square, which can be interpreted with the difference in degrees of freedom of the two models.

.. code-block:: r

        Chi <- LLHom-LLHet
        LRT <- rbind(LLHet,LLHom,Chi)
        LRT

These models may also be specified using paths instead of matrices. See :ref:`multiplegroups-path-specification` for path specification of these models.
