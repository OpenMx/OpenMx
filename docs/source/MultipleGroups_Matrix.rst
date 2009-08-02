Multiple Groups, Matrix Specification
=====================================

An important aspect of structural equation modeling is the use of multiple groups to compare means and covariances structures between any two (or more) data groups, for example males and females, different ethnic groups, ages etc.  Other examples include groups which have different expected covariances matrices as a function of parameters in the model, and need to be evaluated together to estimated together for the parameters to be identified.

Heterogeneity Model
___________________

We will start with a basic example here, building on modeling means and variances in a saturated model.  Assume we have two groups and we want to test whether they have the same mean and covariance structure.  

Data
^^^^

For this example we simulated two datasets ('xy1' and 'xy2') each with zero means and unit variances, one with a correlation of .5, and the other with a correlation of .4 with 1000 subjects each.  See attached R code for simulation and data summary.

.. code-block:: r

    #Simulate Data
    require(MASS)
    #group 1
    set.seed(200)
    rs=.5
    xy1 <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
    set.seed(200)
    #group 2
    rs=.4
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

We first fit a heterogeneity model, allowing differences in both the mean and covariance structure of the two groups.  As we are interested whether the two structures can be equated, we have to specify the models for the two groups, named 'group1' and 'group2' within another model, named 'bivHet'.  The structure of the job thus look as follows, with two ``mxModel`` commands as arguments of another ``mxModel`` command.  ``mxModel`` commands are unlimited in the number of arguments.

.. code-block:: r

    bivHetModel <- mxModel("bivHet",
         mxModel("group1", ....
         mxModel("group2", ....
         mxAlgebra(group1.objective + group2.objective, name="h12"),
         mxAlgebraObjective("h12")
         )
     
For each of the groups, we fit a saturated model, using a Cholesky decomposition to generate the expected covariance matrix and a row vector for the expected means.  Note that we have specified different labels for all the free elements, in the two ``mxModel`` statements.  For more details, see example 1.

.. code-block:: r

    #Fit Heterogeneity Model
    bivHetModel <- mxModel("bivHet",
        mxModel("group1",
            mxMatrix(
                type="Full", 
                nrow=2, 
                ncol=2, 
                free=c(T,T,F,T), 
                values=c(1,.2,0,1),
                labels=c("vX1", "cXY1", "zero", "vY1"),
                dimnames=list(selVars, selVars), 
                name="Chol1"
            ), 
            mxAlgebra(
                Chol1 %*% t(Chol1), 
                name="EC1", 
                dimnames=list(selVars, selVars)
            ), 
            mxMatrix(
                type="Full", 
                nrow=1, 
                ncol=2, 
                free=T, 
                values=c(0,0), 
                labels=c("mX1", "mY1"), 
                dimnames=list(NULL, selVars), 
                name="EM1"
            ), 
            mxData(
                xy1, 
                type="raw"
            ), 
            mxFIMLObjective(
                "EC1", 
                "EM1")
            ),
        mxModel("group2",
            mxMatrix(
                type="Full", 
                nrow=2, 
                ncol=2, 
                free=c(T,T,F,T), 
                values=c(1,.2,0,1),
                labels=c("vX2", "cXY2", "zero", "vY2"),
                dimnames=list(selVars, selVars), 
                name="Chol2"
            ), 
            mxAlgebra(
                Chol2 %*% t(Chol2), 
                name="EC2", 
                dimnames=list(selVars, selVars)
            ), 
            mxMatrix(
                type="Full", 
                nrow=1, 
                ncol=2, 
                free=T, 
                values=c(0,0), 
                labels=c("mX2", "mY2"), 
                dimnames=list(NULL, selVars), 
                name="EM2"
            ), 
            mxData(
                xy2, 
                type="raw"
            ), 
            mxFIMLObjective(
                "EC2", 
                "EM2")
            ), ....

As a result, we estimate five parameters (two means, two variances, one covariance) per group for a total of 10 free parameters.  We cut the 'Labels matrix:' parts from the output generated with ``bivHetModel$group1@matrices`` and ``bivHetModel$group2@matrices``

.. code-block:: r

            $Chol1
              X      Y     
            X "vX1"  "zero"
            Y "cXY1" "vY1" 

            $EM1
                 X     Y    
            [1,] "mX1" "mY1"

            $Chol2
              X      Y     
            X "vX2"  "zero"
            Y "cXY2" "vY2" 

            $EM2
                 X     Y    
            [1,] "mX2" "mY2"

Model Fitting
^^^^^^^^^^^^^

To evaluate both models together, we use an ``mxAlgebra`` command that adds up the values of the objective functions of the two groups.  The objective function to be used here is the ``mxAlgebraObjective`` which uses as its argument the sum of the function values of the two groups.

.. code-block:: r

        mxAlgebra(
                group1.objective + group2.objective, 
                name="h12"
            ),
        mxAlgebraObjective("h12")
        )

The ``mxRun`` command is required to actually evaluate the model.  Note that we have adopted the following notation of the objects.  The result of the ``mxModel`` command ends in 'Model'; the result of the ``mxRun`` command ends in 'Fit'.  Of course, these are just suggested naming conventions.

.. code-block:: r

    bivHetFit <- mxRun(bivHetModel)

A variety of output can be printed.  We chose here to print the expected means and covariance matrices for the two groups and the likelihood of data given the model.  The ``mxEvaluate`` command takes any R expression, followed by the fitted model name.  Given that the model 'bivHetFit' included two models (group1 and group2), we need to use the two level names, i.e. 'group1.EM1' to refer to the objects in the correct model.

.. code-block:: r
    
        EM1Het <- mxEvaluate(group1.EM1, bivHetFit)
        EM2Het <- mxEvaluate(group2.EM2, bivHetFit)
        EC1Het <- mxEvaluate(group1.EC1, bivHetFit)
        EC2Het <- mxEvaluate(group2.EC2, bivHetFit)
        LLHet <- mxEvaluate(objective, bivHetFit)


Homogeneity Model: a Submodel
-----------------------------

Next, we fit a model in which the mean and covariance structure of the two groups are equated to one another, to test whether there are significant differences between the groups.  Rather than having to specify the entire model again, we copy the previous model 'bivHetModel' into a new model 'bivHomModel' to represent homogeneous structures.

.. code-block:: r

    #Fit Homnogeneity Model
    bivHomModel <- bivHetModel

As elements in matrices can be equated by assigning the same label, we now have to equate the labels of the free parameters in group1 to the labels of the corresponding elements in group2.  This can be done by referring to the relevant matrices using the ``ModelName[['MatrixName']]`` syntax, followed by ``@labels``.  Note that in the same way, one can refer to other arguments of the objects in the model.  Here we assign the labels from group1 to the labels of group2, separately for the Cholesky matrices used for the expected covariance matrices and for the expected means vectors.

.. code-block:: r

        bivHomModel[['group2.Chol2']]@labels <- bivHomModel[['group1.Chol1']]@labels
        bivHomModel[['group2.EM2']]@labels <- bivHomModel[['group1.EM1']]@labels

The specification for the submodel is reflected in the names of the labels which are now equal for the corresponding elements of the mean and covariance matrices, as below.

.. code-block:: r

            $Chol1
              X      Y     
            X "vX1"  "zero"
            Y "cXY1" "vY1" 

            $EM1
                 X     Y    
            [1,] "mX1" "mY1"

            $Chol2
              X      Y     
              X "vX1"  "zero"
              Y "cXY1" "vY1" 

            $EM2
                 X     Y    
            [1,] "mX1" "mY1"

We can produce similar output for the submodel, i.e. expected means and covariances and likelihood, the only difference in the code being the model name.  Note that as a result of equating the labels, the expected means and covariances of the two groups should be the same.

.. code-block:: r

    bivHomFit <- mxRun(bivHomModel)
        EM1Hom <- mxEvaluate(group1.EM1, bivHomFit)
        EM2Hom <- mxEvaluate(group2.EM2, bivHomFit)
        EC1Hom <- mxEvaluate(group1.EC1, bivHomFit)
        EC2Hom <- mxEvaluate(group2.EC2, bivHomFit)
        LLHom <- mxEvaluate(objective, bivHomFit)

Finally, to evaluate which model fits the data best, we generate a likelihood ratio test as the difference between -2 times the log-likelihood of the homogeneity model and -2 times the log-likelihood of the heterogeneity model.  This statistic is asymptotically distributed as a Chi-square, which can be interpreted with the difference in degrees of freedom of the two models.

.. code-block:: r

        Chi= LLHom-LLHet
        LRT= rbind(LLHet,LLHom,Chi)
        LRT
