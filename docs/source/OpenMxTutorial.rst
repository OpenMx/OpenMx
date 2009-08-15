Tutorial
========

Prerequisites
-------------

Congratulations!  You have decided to check out OpenMx, the open source version of the statistical modeling package Mx, rewritten in R.  Before we get started, let's make sure you have the software installed and ready to go.  As
You need:

* `R <http://www.r-project.org/>`_
* OpenMx_

.. _OpenMx: http://openmx.psyc.virginia.edu


Simple OpenMx Script
--------------------

We will start by showing some of the main features of OpenMx using simple examples.  For those familiar with Mx, it is basically a matrix interpreter combined with a numerical optimizer to allow fitting statistical models.  Of course you do not need OpenMx to perform matrix algebra as that can already be done in R.  However, to accomodate flexible statistical modeling of the type of models typically fit in Mx, Mplus or other SEM packages, special kinds of matrices and functions are required which are bundled in OpenMx.  We will introduce key features of OpenMx using a matrix algebra example.  Remember that R is object-oriented, such that the results of operations are objects, rather than just matrices, with various properties/characteristics attached to them.

Say, we want to create two matrices, **A** and **B**, each of them a 'Full' matrix with 3 rows and 1 column and with the values 1, 2 and 3, as follows:

.. math::
   :nowrap:

   \begin{eqnarray*}
   A = \left[ \begin{array}{r} 1 \\ 2 \\ 3 \\ \end{array} \right]
   & B = \left[ \begin{array}{r} 1 \\ 2 \\ 3 \\ \end{array} \right]
   \end{eqnarray*}

we use the ``mxMatrix`` command, and define the type of the matrix, number of rows and columns, its specifications and values, <optionally labels, upper and lower bounds>,  and a name.  The matrix **A** will be stored as the object 'A'.

.. code-block:: r

    mxMatrix(
        type="Full", 
        nrow=3, 
        ncol=1, 
        values=c(1,2,3), 
        name='A'
    )
    mxMatrix(
        type="Full", 
        nrow=3, 
        ncol=1, 
        values=c(1,2,3), 
        name='A'
    )

Assume we want to calculate	the (1) outer and (2) inner products of the matrix **A**, using regular matrix multiplication, an (3) element by element multiplication (Dot product) and (4) the sum of the matrices **A** and **B**, i.e.:

.. math::
   :nowrap:

    \begin{eqnarray}
    q1 & = & A * t(A) \\
    q2 & = & t(A) * A \\
    q3 & = & A . A \\
    q4 & = & A + B
    \end{eqnarray}

we invoke the ``mxAlgebra`` command which takes an algebra operation between previously defined matrices.  Note that in R, regular matrix multiplication is represented by ``\%*\%`` and dot multiplication as ``*``. We also assign the algebras a name to refer back to them later:

.. code-block:: r

    mxAlgebra(
        A %*% t(A), 
        name='q1'
    )
    mxAlgebra(
        t(A) %*% A, 
        name='q2'
    )
    mxAlgebra(
        A * A, 
        name='q3'
    )
    mxAlgebra(
        A + B, 
        name='q4'
    )

For the algebras to be evaluated, they become arguments of the ``mxModel`` command, as do the defined matrices, separated by comma's.  The model, which is here given the name 'algebraExercises', is then executed by the ``mxRun`` command, as shown in the full code below:

.. code-block:: r

    require(OpenMx)

    algebraExercises <- mxModel(
        mxMatrix(type="Full", values=c(1,2,3), nrow=3, ncol=1, name='A'),
        mxMatrix(type="Full", values=c(1,2,3), nrow=3, ncol=1, name='B'),
        mxAlgebra(A%*%t(A), name='q1'),
        mxAlgebra(t(A)%*%A, name='q2'),
        mxAlgebra(A*A, name='q3'),
        mxAlgebra(A+B, name='q4'))

    answers <- mxRun(algebraExercises)
    answers@algebras
    result <- mxEval(list(q1,q2,q3,q4),answers)	

As you notice, we added some lines at the end to generate the desired output.  The resulting matrices and algebras are stored in ``answers``; we can refer back to them by specifying ``answers@matrices`` or ``answers@algebras``.  We can also calculate any additional quantities or perform extra matrix operations on the results using the ``mxEval`` command.  For example, if we want to see all the answers to the questions in matrixAlgebra.R, the results would look like this:

.. code-block:: r

    [[1]]
         [,1] [,2] [,3]
    [1,]    1    2    3
    [2,]    2    4    6
    [3,]    3    6    9

    [[2]]
         [,1]
    [1,]   14

    [[3]]
         [,1]
    [1,]    1
    [2,]    4
    [3,]    9

    [[4]]
         [,1]
    [1,]    2
    [2,]    4
    [3,]    6

So far, we have introduced five new commands: ``mxMatrix``, ``mxAlgebra``, ``mxModel``, ``mxRun`` and ``mxEval``.  These commands allow us to run a wide range of jobs, from simple matrix algebra to rather complicated SEM models.  Let's move to a simple example involving optimizing the likelihood of observed data.


Optimization Script
-------------------

When collecting data to test a specific hypothesis, one of the first things one typically does is checking the basic descriptive statistics, such as the means, variances and covariances/correlations.  We could of course use basic functions in R, i.e., `meanCol(Data)` or `cov(Data)`.  However, if we want to test specific hypotheses about the data, for example, test whether the correlation between our two variables is significantly different from zero, we need to compare the likelihood of the data with that where the correlation is fixed to zero.  Let's work through a specific example.

Say, we have collected data on two variables **X** and **Y** in 1000 individuals, and R descriptive statistics has shown that the correlation between them in 0.5.  For the sake of this example, we used another built-in function in the R package MASS, namely mvrnorm, to generate multivariate normal data with means of 0.0, variances of 1.0 and a correlation of 0.5 between **X** and **Y**.

To evaluate the likelihood of the data, we estimate a saturated model with free means, free variances and a covariance.  Let's start with specifying the mean vector.  We use the ``mxMatrix`` command, provide the type, here "Full", the number of rows and columns, respectively 1 and 2, the specification of free/fixed parameters, the starting values, the dimnames and a name.  Given all the elements of this 1x2 matrix are free, we can use ``free=True``.  The starting values are provided using a list, i.e. ``c(0,0)``.  The dimnames are a type of label that is required to recognize the expected mean vector and expected covariance matrix.  In this case, the second element of the list should have the labels for the two variables ``c('X','Y')``.  Finally, we are explicit in naming this matrix ``expMean``.  Thus the matrix command looks like this.  Note the soft tabs to improve readability.

.. code-block:: r

    bivCorModel <- mxModel("bivCor",
        mxMatrix(
            type="Full", 
            nrow=1, 
            ncol=2, 
            free=True, 
            values=c(0,0), 
            dimnames=list(NULL, selVars), 
            name="expMean"
        ), 

Next, we need to specify the expected covariance matrix.  As this matrix is symmetric, we could estimate it directly as a symmetric matrix.  However, to avoid solutions that are not positive definite, we will use a Cholesky decomposition.  Thus, we specify a lower triangular matrix (matrix with free elements on the diagonal and below the diagonal, and zero's above the diagonal), and multiply it with its transpose to generate a symmetric matrix.  We will use a ``mxMatrix`` command to specify the lower triangular matrix and a ``mxAlgebra`` command to set up the symmetric matrix.  (PS a lower triangular matrix doesn't exist yet so we specify it explicitly.)  The matrix is a 2x2 free lower matrix with  ``c('X','Y')`` as dimnames for the rows and columns, and the name "Chol".  We can now refer back to this matrix by its name in the ``mxAlgebra`` statement.  We use a regular multiplication of ``Chol`` with its transpose ``t(Chol)``, and name this as "expCov".

.. code-block:: r

        mxMatrix(
            type="Full", 
            nrow=2, 
            ncol=2, 
            free=c(T,T,F,T), 
            values=c(1,.2,0,1), 
            dimnames=list(selVars, selVars), 
            name="Chol"
        ), 
        mxAlgebra(
            expression=Chol %*% t(Chol), 
            name="expCov", 
            dimnames=list(selVars, selVars)
        ), 

Now that we have specified our 'model', we need to supply the data.  This is done with the ``mxData`` command.  The first argument includes the actual data, in the type given by the second argument.  Type can be a covariance matrix (cov), a correlation matrix (cor), a matrix of cross-products (sscp) or raw data (raw).  We will use the latter option and read in the raw data directly from the simulated dataset ``testData``.

.. code-block:: r

        mxData(
            observed=testData, 
            type="raw"
        ), 

Next, we specify which objective function we wish to use to obtain the likelihood of the data.  Given we fit to the raw data, we use the full information maximum likelihood (FIML) objective function ``mxFIMLObjective``.  Its arguments are the expected covariance matrix, generated using the ``mxMatrix`` and ``mxAlgebra`` commands as "expCov", and the expected means vectors, generated using the ``mxMatrix`` command as "expMeans".

.. code-block:: r

        mxFIMLObjective(
            covariance="expCov", 
            means="expMean")
        )

All these elements become arguments of the ``mxModel`` command, seperated by comma's.  The first argument can be a name, as in this case "bivCor" or another model (see below).  The model is then saved in an object 'bivCorModel' which becomes the argument of the ``mxRun`` command, which evaluates the model and provides output - if the model ran successfully. using the following command.

.. code-block:: r

        bivCorFit <- mxRun(bivCorModel)

We can then request various parts of the output to inspect by referring to them by the name of the object resulting from the ``mxRun`` command, followed by the name of the objects corresponding to the expected mean vector and covariance matrix, in quotes and double square brackets, followed by ``@values``.  The command ``mxEval`` can also be used to extract relevant information, such as the likelihood, where the first argument of the command is the object of interest and the second the object obtaining the results.

.. code-block:: r

    EM <- bivCorFit[['expMean']]@values
    EC <- bivCorFit[['expCov']]@values
    LL <- mxEval(objective,bivCorFit);

If we want to test whether the covariance/correlation is significantly different from zero, we could fit a submodel and compare it with the saturated model.  Given that this model is essentially the same as the original, except for the covariance, we create a new mxModel (named 'bivCorModelSub) with as first argument the old model (named 'bivCorModel).  Then we only have to specify the matrix that needs to be changed, in this case the lower triangular matrix becomes essentially a diagonal matrix, obtained by fixing the off-diagonal elements to zero in the ``free`` and ``values`` arguments

.. code-block:: r

    #Test for Covariance=Zero
    bivCorModelSub <-mxModel(bivCorModel,
        mxMatrix(
            type="Full", 
            nrow=2, 
            ncol=2, 
            free=c(T,F,F,T), 
            values=c(1,0,0,1), 
            dimnames=list(selVars, selVars),
            name="Chol"
        )

We can output the same information as for the saturated job, namely the expected means and covariance matrix and the likelihood, and then use R to calculate other statistics, such as the Chi-square goodness-of-fit.

.. code-block:: r

    bivCorFitSub <- mxRun(bivCorModelSub)
    EMs <- bivCorFitSub[['expMean']]@values
    ECs <- bivCorFitSub[['expCov']]@values
    LLs <- mxEval(objective,bivCorFitSub);
    Chi= LLs-LL;
    LRT= rbind(LL,LLs,Chi); LRT


More in-depth Example
---------------------

Now that you have seen the basics of OpenMx, let us walk through an example in more detail.  We decided to use a twin model example for several reasons.  Even though you may not have any background in behavior genetics or genetic epidemiology, the example illustrates a number of features you are likely to encounter at some stage.  We will present the example in two ways: (i) path analysis representation, and (ii) matrix algebra representation.  Both give exactly the same answer, so you can choose either one or both to get some familiarity with the two approaches.

We will not go into detail about the theory of this model, as that has been done elsewhere (refs).  Briefly, twin studies rely on comparing the similarity of identical (monozygotic, MZ) and fraternal (dizygotic, DZ) twins to infer the role of genetic and environmental factors on individual differences.  As MZ twins have identical genotypes, similarity between MZ twins is a function of shared genes, and shared environmental factors.  Similarity between DZ twins is a function of some shared genes (on average they share 50% of their genes) and shared environmental factors.  A basic assumption of the classical twin design is that the MZ and DZ twins shared environmental factors to the same extent.

The basic model typically fit to twin data from MZ and DZ twins reared together includes three sources of latent variables: additive genetic factors (**A**), shared environmental influences (**C**) and unique environmental factors (**E**),  We can estimate these three sources of variance from the observed variances, the MZ and the DZ covariance.  The expected variance is the sum of the three variance components (**A + C + E**).  The expected covariance for MZ twins is (**A + C**) and that of DZ twins is (**.5A + C**).  As MZ and DZ twins have different expected covariances, we have multiple group model.

It has been standard in twin modeling to fit models to the raw data, as often data are missing on some co-twins.  When using FIML, we also need to specify the expected means.  There is no reason to expect that the variances are different for twin 1 and twin 2, neither are the means for twin 1 and twin 2 expected to differ.  This can easily be verified by fitting submodels to the saturated model, prior to fitting the ***ACE*** model.

Let us start by fitting a saturated model, estimating means, variances and covariances separately order of the twins (twin 1 vs twin 2) and by zygosity (MZ vs DZ pairs).  This is essentially similar to the optimization script discussed above, except that we now have two variables (same variable for twin 1 and twin 2) and two groups (MZ and DZ).  Before we get to the OpenMx code, let us organize the data in R.

.. code-block:: r

    require(OpenMx)

    #Prepare Data
    twinData <- read.table("myTwinData.txt", header=T, na.strings=".")
    twinVars <- c('fam','age','zyg','part','wt1','wt2','ht1','ht2','htwt1','htwt2','bmi1','bmi2')
    summary(twinData)
    selVars <- c('bmi1','bmi2')
    mzfData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2)))
    dzfData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2)))

The saturated model will have two matrices for the expected means of MZs and DZs, and two for the expected covariances, generated from multiplying a lower triangular matrix with its transpose.  The raw data are read in using the ``mxData`` command, and the corresponding objective funtion ``mxFIMLObjective`` applied.  

.. code-block:: r

    mxModel("MZ",
        mxMatrix(
            type="Full", 
            nrow=1, 
            ncol=2, 
            free=True, 
            values=c(0,0), 
            dimnames=list(NULL, selVars), 
            name="expMeanMZ"), 
        mxMatrix("Full", 2, 2,
            free=c(T,T,F,T)
            values=c(1,.5,0,1), 
            dimnames=list(NULL, selVars), 
            name="CholDZ"), 
        mxAlgebra(
            CholMZ %*% t(CholMZ), 
            name="expCovMZ", 
            dimnames=list(selVars, selVars)), 
        mxData(
            DataMZ, 
            type="raw"), 
        mxFIMLObjective(
            "expCovMZ", 
            "expMeanMZ"))

Note that the ``mxModel`` statement for the DZ twins is almost identical to that for MZ twins, except for the names of the objects and data.  If the arguments to the OpenMx command are given in the default order (see i.e. ?mxMatrix to go to the help/reference page for that command), then it is not necessary to include the name of the argument.  Given we skip a few optional arguments, ``dimnames=`` and ``name=`` are included to refer to the right arguments.  For didactic purposes, we prefer the formatting used for the MZ group, with soft tabs and each argument on a separate line, etc.  (see list of formatting rules).  However, the experienced user may want to use a more compact form, as the one used for the DZ group.

.. code-block:: r            

    mxModel("DZ",
        mxMatrix("Full", 1, 2, T, c(0,0), dimnames=list(NULL, selVars), name="expMeanDZ"), 
        mxMatrix("Full", 2, 2, c(T,T,F,T), c(1,.5,0,1), dimnames=list(NULL, selVars), name="CholMZ"), 
        mxAlgebra(CholDZ %*% t(CholDZ), name="expCovDZ", dimnames=list(selVars, selVars)), 
        mxData(DataDZ, type="raw"), 
        mxFIMLObjective("expCovDZ", "expMeanDZ")),

The two models are then combined in a 'super'model which includes them as arguments.  Additional arguments are an ``mxAlgebra`` statement to add the objective funtions/likelihood of the two submodels.  To evaluate them simultaneously, we use the ``mxAlgebraObjective`` with the previous algebra as its argument.  The ``mxRun`` command is used to start optimization.

.. code-block:: r 

    twinSatModel <- mxModel("twinSat",
        mxModel("MZ", .... ),
        mxModel("DZ", .... ),
        mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
        mxAlgebraObjective("twin"))
    twinSatFit <- mxModel(twinSatModel)

It is always helpful/advised to check the model specifications before interpreting the output.  Here we are interested in the values for the expected mean vectors and covariance matrices, and the goodness-of-fit statistics, including the likelihood, degrees of freedom, and any other derived indices.

.. code-block:: r

    ExpMeanMZ <- mxEval(MZ.expMeanMZ, twinSatFit)
    ExpCovMZ <- mxEval(MZ.expCovMZ, twinSatFit)
    ExpMeanDZ <- mxEval(DZ.expMeanDZ, twinSatFit)
    ExpCovDZ <- mxEval(DZ.expCovDZ, twinSatFit)
    LL_Sat <- mxEval(objective, twinSatFit)

Before we move on to fit the ACE model to the same data, we may want to test some of the assumptions of the twin model, i.e. that the means and variances are the same for twin 1 and twin 2, and that they are the same for MZ and DZ twins.  This can be done as an omnibus test, or stepwise.  Let us start by equating the means for both twins, separately in the two groups.  As the majority of the previous script stays the same, we start by copying the old model into a new one.  We then include the arguments of the model that require a change.

.. code-block:: r 

    twinSatModelSub1 <- mxModel(twinSatModel,
        mxModel("MZ",
            mxMatrix("Full", 1, 2, T, 0, "mMZ", dimnames=list(NULL, selVars), name="expMeanMZ"), 
        mxModel("DZ", 
            mxMatrix("Full", 1, 2, T, 0, "mDZ", dimnames=list(NULL, selVars), name="expMeanDZ"), 
        mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
        mxAlgebraObjective("twin"))
    twinSatFitSub1 <- mxModel(twinSatModelSub1)

If we want to test if we can equate both means and variances across twin order and zygosity at once, we will end up with the following specification.  Note that we use the same label for elements that need to be equated.

.. code-block:: r 

    twinSatModelSub2 <- mxModel(twinSatModelSub1,
        mxModel("MZ",
            mxMatrix("Full", 1, 2, T, 0, "mean", dimnames=list(NULL, selVars), name="expMeanMZ"), 
            mxMatrix("Full", 2, 2, c(T,T,F,T), c(1,.5,0,1), labels= c("var","MZcov","var"), 
                dimnames=list(NULL, selVars), name="CholMZ"), 
        mxModel("DZ", 
            mxMatrix("Full", 1, 2, T, 0, "mean", dimnames=list(NULL, selVars), name="expMeanDZ"), 
            mxMatrix("Full", 2, 2, c(T,T,F,T), c(1,.5,0,1), labels= c("var","DZcov","var"), 
                dimnames=list(NULL, selVars), name="CholDZ"), 
        mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
        mxAlgebraObjective("twin"))
    twinSatFitSub2 <- mxModel(twinSatModelSub2)

We can compare the likelihood of this submodel to that of the fully saturated model or the previous submodel using the results from ``mxEval`` commands with regular R algebra.  A summary of the model parameters, estimates and goodness-of-fit statistics can also be obtained using ``summary(twinSatFit)``.  Further development is required.

.. code-block:: r

    LL_Sat <- mxEval(objective, twinSatFit)
    LL_Sub <- mxEval(objective, twinSatFitSub1);
    LRT= LL_Sub - LL_Sat;

Now, we are ready to specify the ACE model to test which sources of variance significantly contribute to the phenotype and estimate their best value.  The structure of this script is going to mimic that of the saturated model.  The main difference is that we no longer estimate the variance-covariance matrix directly, but express it as a function of the three sources of variance, **A**, **C** and **E**.  As the same sources are used for the MZ and the DZ group, the matrices which will represent them are part of the 'super'model.  As these sources are variances, which need to be positive, we typically use a Cholesky decomposition of the standard deviations (and effectively estimate **a** rather then **a^2**, see later for more in depth coverage).  Thus, we specify three separate matrices for the three sources of variance using the ``mxMatrix`` command and 'calculate' the variance components with the ``mxAlgebra`` command.  Note that there are a variety of ways to specify this model, we have picked one that corresponds well to previous Mx code, and has some intuitive appeal.

.. code-block:: r

    #Specify ACE Model
    twinACEModel <- mxModel("twinACE", 
        mxMatrix("Full", 1, 2, T, 20, "mean", dimnames=list(NULL, selVars), name="expMeanMZ"), 
        mxMatrix("Full", 1, 2, T, 20, "mean", dimnames=list(NULL, selVars), name="expMeanDZ"), 
        mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="a", name="X"),
        mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="c", name="Y"),
        mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z"),
        mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.5, name="h"),
        mxAlgebra(X * t(X), name="A"),
        mxAlgebra(Y * t(Y), name="C"),
        mxAlgebra(Z * t(Z), name="E"), 
        mxAlgebra(rbind (cbind(A+C+E   , A+C),
                         cbind(A+C     , A+C+E)), dimnames = list(selVars, selVars), name="expCovMZ"),
        mxAlgebra(rbind (cbind(A+C+E   , h%x%A+C),
                         cbind(h%x%A+C , A+C+E)), dimnames = list(selVars, selVars), name="expCovDZ"),
        mxModel("MZ",
            mxData(mzfData, type="raw"), 
            mxFIMLObjective("twinACE.expCovMZ", "twinACE.expMeanMZ")),
        mxModel("DZ", 
            mxData(dzfData, type="raw"), 
            mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMeanDZ")),
        mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
        mxAlgebraObjective("twin"))
    #Run ACE model can be run
    twinACEFit <- mxRun(twinACEModel)

Relevant output can be generate with ``print`` or ``summary`` statements or specific output can be requested using the ``mxEval`` command.  Typically we would compare this model back to the saturated model to interpret its goodness-of-fit.  Parameter estimates are obtained and can easily be standardized.  We discuss a twin analysis example in more detail in the example code.