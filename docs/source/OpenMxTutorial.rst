Tutorial
========

Prerequisites
-------------

Congratulations!  You have decided to check out OpenMx, the open source version of the statistical modeling package Mx, rewritten in R.  Before we get started, let's make sure you have the software installed and ready to go.  
You need:

* `R <http://www.r-project.org/>`_
* OpenMx_

.. _OpenMx: http://openmx.psyc.virginia.edu

------------


Simple OpenMx Script
--------------------

We will start by showing some of the main features of OpenMx using simple examples.  For those familiar with Mx, it is basically a matrix interpreter combined with a numerical optimizer to allow fitting statistical models.  Of course you do not need OpenMx to perform matrix algebra as that can already be done in R.  However, to accommodate flexible statistical modeling of the type of models typically fit in Mx, Mplus or other SEM packages, special kinds of matrices and functions are required which are bundled in OpenMx.  We will introduce key features of OpenMx using a matrix algebra example.  Remember that R is object-oriented, such that the results of operations are objects, rather than just matrices, with various properties/characteristics attached to them.  We will describe the script line by line; a link to the complete script is `here
<http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/MatrixAlgebra20090924.R>`_.

Say, we want to create two matrices, **A** and **B**, each of them a 'Full' matrix with 3 rows and 1 column and with the values 1, 2 and 3, as follows:

.. math::
   :nowrap:

   \begin{eqnarray*}
   A = \left[ \begin{array}{r} 1 \\ 2 \\ 3 \\ \end{array} \right]
   & B = \left[ \begin{array}{r} 1 \\ 2 \\ 3 \\ \end{array} \right]
   \end{eqnarray*}

we use the ``mxMatrix`` command, and define the type of the matrix (``type=``), number of rows (``nrow=``) and columns (``ncol=``), its specifications (``free=``) and starting values (``values=``), optionally labels (``labels=``), upper (``ubound=``) and lower (``lbound=``) bounds>, and a name (``name=``).  The matrix **A** will be stored as the object 'A'.

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
        name='B'
    )

Assume we want to calculate	the (1) the sum of the matrices **A** and **B**, (2) the element by element multiplication (Dot product) of **A** and **B**, (3) the transpose of matrix **A**, and the (4) outer and (5) inner products of the matrix **A**, using regular matrix multiplication, i.e.:

.. math::
   :nowrap:

    \begin{eqnarray}
    q2 & = & A + B \\
    q1 & = & A . A \\
	q3 & = & t(A) \\
    q4 & = & A * t(A) \\
    q5 & = & t(A) * A
    \end{eqnarray}

we invoke the ``mxAlgebra`` command which performs an algebra operation between previously defined matrices.  Note that in R, regular matrix multiplication is represented by ``\%*\%`` and dot multiplication as ``*``. We also assign the algebras a name to refer back to them later:

.. code-block:: r

	mxAlgebra(
		A + B, 
		name='q1'
	)
	mxAlgebra(
		A * A, 
		name='q2'
	)
	mxAlgebra(
		t(A), 
		name='q3'
	)
	mxAlgebra(
		A %*% t(A), 
		name='q4'
	)
	mxAlgebra(
		t(A) %*% A, 
		name='q5'
	)

For the algebras to be evaluated, they become arguments of the ``mxModel`` command, as do the defined matrices, each separated by comma's.  The model, which is here given the name 'algebraExercises', is then executed by the ``mxRun`` command, as shown in the full code below:

.. code-block:: r

    require(OpenMx)

    algebraExercises <- mxModel(
        mxMatrix(type="Full", nrow=3, ncol=1, values=c(1,2,3), name='A'),
        mxMatrix(type="Full", nrow=3, ncol=1, values=c(1,2,3), name='B'),      
        mxAlgebra(A+B, name='q1'),
		mxAlgebra(A*A, name='q2'),
		mxAlgebra(t(A), name='q3'),		
        mxAlgebra(A%*%t(A), name='q4'),
        mxAlgebra(t(A)%*%A, name='q5'))

    answers <- mxRun(algebraExercises)
    answers@algebras
    result <- mxEval(list(q1,q2,q3,q4,q5),answers)	

As you notice, we added some lines at the end to generate the desired output.  The resulting matrices and algebras are stored in ``answers``; we can refer back to them by specifying ``answers@matrices`` or ``answers@algebras``.  We can also calculate any additional quantities or perform extra matrix operations on the results using the ``mxEval`` command.  For example, if we want to see all the answers to the questions in matrixAlgebra.R, the results would look like this:

.. code-block:: r

	[[1]]
		[,1]
	[1,]    2
	[2,]    4
	[3,]    6

	[[2]]
		[,1]
	[1,]    1
	[2,]    4
	[3,]    9

	[[3]]
	  	[,1] [,2] [,3]
	[1,]    1    2    3

	[[4]]
		[,1] [,2] [,3]
	[1,]    1    2    3
	[2,]    2    4    6
	[3,]    3    6    9

	[[5]]
		[,1]
	[1,]   14


So far, we have introduced five new commands: ``mxMatrix``, ``mxAlgebra``, ``mxModel``, ``mxRun`` and ``mxEval``.  These commands allow us to run a wide range of jobs, from simple matrix algebra to rather complicated SEM models.  Let's move to an example involving optimizing the likelihood of observed data.

------------

Optimization Script
-------------------

When collecting data to test a specific hypothesis, one of the first things one typically does is checking the basic descriptive statistics, such as the means, variances and covariances/correlations.  We could of course use basic functions in R, i.e., `meanCol(Data)` or `cov(Data)` to perform these operations.  However, if we want to test specific hypotheses about the data, for example, test whether the correlation between two variables is significantly different from zero, we need to compare the likelihood of the data when the correlation is freely estimated with the likelihood of the data when the correlation is fixed to zero.  Let's work through a specific `example
<http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/BivariateCorrelation20090925.R`_>.

Say, we have collected data on two variables **X** and **Y** in 1000 individuals, and R descriptive statistics has shown that the correlation between them is 0.5.  For the sake of this example, we used another built-in function in the R package MASS, namely ``mvrnorm``, to generate multivariate normal data for 1000 individuals with means of 0.0, variances of 1.0 and a correlation (``rs``) of 0.5 between **X** and **Y**.  Note the that first argument of ``mvrnorm`` is the sample size, the second the vector of means, and the third the covariance matrix to be simulated.  We save the data in the object ``xy`` and create a vector of labels for the two variables in ``selVars`` which is used in the ``dimnames`` statement later on.  The R functions ``summary()`` and ``cov()`` are used to verify that the simulations appear OK.

.. code-block:: r

	#Simulate Data
	require(MASS)
	set.seed(200)
	rs=.5
	xy <- mvrnorm (1000, c(0,0), matrix(c(1,rs,rs,1),2,2))
	testData <- xy
	selVars <- c('X','Y')
	dimnames(testData) <- list(NULL, selVars)
	summary(testData)
	cov(testData)

To evaluate the likelihood of the data using SEM, we estimate a saturated model with free means, free variances and a covariance.  Let's start with specifying the mean vector.  We use the ``mxMatrix`` command, provide the ``type``, here ``Full``, the number of rows and columns, respectively 1 and 2, the specification of free/fixed parameters, the starting values, the dimnames and a name.  Given all the elements of this 1x2 matrix are free, we can use ``free=True``.  The starting values are provided using a list, i.e. ``c(0,0)``.  The ``dimnames`` are a type of label that is required to recognize the expected mean vector and expected covariance matrix and match up the model with the data.  For a mean vector, the first element is ``NULL`` given mean vectors always have one row.  The second element of the list should have the labels for the two variables ``c('X','Y')`` which we have previously assigned to the object ``selVars``.  Finally, we are explicit in naming this matrix ``expMean``.  Thus the matrix command looks like this.  Note the soft tabs to improve readability.

.. code-block:: r

    bivCorModel <- mxModel("bivCor",
        mxMatrix(
            type="Full", 
            nrow=1, 
            ncol=2, 
            free=TRUE, 
            values=c(0,0), 
            name="expMean"
        ), 

Next, we need to specify the expected covariance matrix.  As this matrix is symmetric, we could estimate it directly as a symmetric matrix.  However, to avoid solutions that are not positive definite, we will use a Cholesky decomposition.  Thus, we specify a lower triangular matrix (matrix with free elements on the diagonal and below the diagonal, and zero's above the diagonal), and multiply it with its transpose to generate a symmetric matrix.  We will use a ``mxMatrix`` command to specify the lower triangular matrix and a ``mxAlgebra`` command to set up the symmetric matrix.  The matrix is a 2x2 free lower matrix with ``c('X','Y')`` (previously defined as ``selVars``) as ``dimnames`` for the rows and columns, and the name "Chol".  We can now refer back to this matrix by its name in the ``mxAlgebra`` statement.  We use a regular multiplication of ``Chol`` with its transpose ``t(Chol)``, and name this as "expCov".

.. code-block:: r

        mxMatrix(
            type="Lower", 
            nrow=2, 
            ncol=2, 
            free=TRUE, 
            values=.5, 
            name="Chol"
        ), 
        mxAlgebra(
            expression=Chol %*% t(Chol), 
            name="expCov"
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
            means="expMean",
            dimnames=selVars)
        )

All these elements become arguments of the ``mxModel`` command, seperated by comma's.  The first argument can be a name, as in this case "bivCor" or another model (see below).  The model is saved in an object 'bivCorModel'.  This object becomes the argument of the ``mxRun`` command, which evaluates the model and provides output - if the model ran successfully - using the following command.

.. code-block:: r

        bivCorFit <- mxRun(bivCorModel)

We can request various parts of the output to inspect by referring to them by the name of the object resulting from the ``mxRun`` command, i.e. ``bivCorFit``, followed by the name of the objects corresponding to the expected mean vector, i.e. ``[['ExpMean']]``, and covariance matrix, i.e. ``[['ExpCov']]``, in quotes and double square brackets, followed by ``@values``.  The command ``mxEval`` can also be used to extract relevant information, such as the likelihood, (``objective``) where the first argument of the command is the object of interest and the second the object obtaining the results.

.. code-block:: r

    EM <- bivCorFit[['expMean']]@values
    EC <- bivCorFit[['expCov']]@values
    LL <- mxEval(objective,bivCorFit);

These commands generate the following output:

.. code-block:: r

 	EM
               	  X            Y
	[1,] 0.03211646 -0.004883803

	EC
          	  X         Y
	X 1.0092847 0.4813501
	Y 0.4813501 0.9935387

	LL
           	 [,1]
	[1,] 5415.772
	
Standard lists of parameter estimates and goodness-of-fit statistics can also be obtained with the ``summary`` command.

.. code-block:: r

	> summary(bivCorFit)
	       X                   Y            
	 Min.   :-2.942561   Min.   :-3.296159  
	 1st Qu.:-0.633711   1st Qu.:-0.596177  
	 Median :-0.004139   Median :-0.010538  
	 Mean   : 0.032116   Mean   :-0.004884  
	 3rd Qu.: 0.739236   3rd Qu.: 0.598326  
	 Max.   : 4.173841   Max.   : 4.006771  

	  name  matrix row col parameter estimate error estimate
	1 <NA> expMean   1   1        0.032116456     0.02228409
	2 <NA> expMean   1   2       -0.004883803     0.02235021
	3 <NA>    Chol   1   1        1.004631642     0.01575904
	4 <NA>    Chol   2   1        0.479130899     0.02099642
	5 <NA>    Chol   2   2        0.874055066     0.01376876

	Observed statistics:  2000 
	Estimated parameters:  5 
	Degrees of freedom:  1995 
	-2 log likelihood:  5415.772 
	Saturated -2 log likelihood:  
	Chi-Square:   
	p:   
	AIC (Mx):  1425.772 
	BIC (Mx):  -4182.6 
	adjusted BIC: 
	RMSEA:  0

If we want to test whether the covariance/correlation is significantly different from zero, we could fit a submodel and compare it with the previous saturated model.  Given that this model is essentially the same as the original, except for the covariance, we create a new mxModel (named 'bivCorModelSub) with as first argument the old model (named 'bivCorModel).  Then we only have to specify the matrix that needs to be changed, in this case the lower triangular matrix becomes essentially a diagonal matrix, obtained by fixing the off-diagonal elements to zero in the ``free`` and ``values`` arguments

.. code-block:: r

    #Test for Covariance=Zero
    bivCorModelSub <-mxModel(bivCorModel,
        mxMatrix(
			type="Diag", 
			nrow=2, 
			ncol=2,
			free=TRUE,
			dimnames=list(selVars, selVars),
			name="Chol"
        )

We can output the same information as for the saturated job, namely the expected means and covariance matrix and the likelihood, and then use R to calculate other statistics, such as the Chi-square goodness-of-fit.

.. code-block:: r

    bivCorFitSub <- mxRun(bivCorModelSub)
    EMs <- mxEval(expMean, bivCorFitSub)
    ECs <- mxEval(expCov, bivCorFitSub)
    LLs <- mxEval(objective, bivCorFitSub)
    Chi= LLs-LL;
    LRT= rbind(LL,LLs,Chi); LRT

------------

More in-depth Example
---------------------

Now that you have seen the basics of OpenMx, let us walk through an example in more detail.  We decided to use a twin model example for several reasons.  Even though you may not have any background in behavior genetics or genetic epidemiology, the example illustrates a number of features you are likely to encounter at some stage.  We will present the example in two ways: (i) path analysis representation, and (ii) matrix algebra representation.  Both give exactly the same answer, so you can choose either one or both to get some familiarity with the two approaches.

We will not go into detail about the theory of this model, as that has been done elsewhere (refs).  Briefly, twin studies rely on comparing the similarity of identical (monozygotic, MZ) and fraternal (dizygotic, DZ) twins to infer the role of genetic and environmental factors on individual differences.  As MZ twins have identical genotypes, similarity between MZ twins is a function of shared genes, and shared environmental factors.  Similarity between DZ twins is a function of some shared genes (on average they share 50% of their genes) and shared environmental factors.  A basic assumption of the classical twin design is that the MZ and DZ twins shared environmental factors to the same extent.

The basic model typically fit to twin data from MZ and DZ twins reared together includes three sources of latent variables: additive genetic factors (**A**), shared environmental influences (**C**) and unique environmental factors (**E**),  We can estimate these three sources of variance from the observed variances, the MZ and the DZ covariance.  The expected variance is the sum of the three variance components (**A + C + E**).  The expected covariance for MZ twins is (**A + C**) and that of DZ twins is (**.5A + C**).  As MZ and DZ twins have different expected covariances, we have a multiple group model.

It has been standard in twin modeling to fit models to the raw data, as often data are missing on some co-twins.  When using FIML, we also need to specify the expected means.  There is no reason to expect that the variances are different for twin 1 and twin 2, neither are the means for twin 1 and twin 2 expected to differ.  This can easily be verified by fitting submodels to the saturated model, prior to fitting the ***ACE*** model.

Let us start by simulating the data following by fitting a series of models.  The `code
<http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/UnivariateTwinAnalysis20090925.R>`_.
includes both the twin data simulation and several OpenMx scripts to analyze the data.  We will describe each of the parts in turn and include the code for the specific part in the code blocks.

First, we simulate twin data using the ``mvrnorm`` R function.  If the additive genetic factors (**A**) account for 50% of the total variance and the shared environmental factors (**C**) for 30%, thus leaving 20% explained by specific environmental factors (**E**), then the expected MZ twin correlation is ``a^2 + c^2`` or 0.8 in this case, and the expected DZ twin correlation is 0.65, calculated as ``.5*a^2 + c^2``.  We simulate 1000 pairs of MZ and DZ twins each with zero means and a correlation matrix according to the values listed above.  We run some basic descriptive statistics on the simulated data, using regular R functions.

.. code-block:: r

    require(OpenMx)

 	require(MASS)
	set.seed(200)
	a2<-0.5		#Additive genetic variance component (a squared)
	c2<-0.3		#Common environment variance component (c squared)
	e2<-0.2		#Specific environment variance component (e squared)
	rMZ <- a2+c2
	rDZ <- .5*a2+c2
	MZ <- mvrnorm (1000, c(0,0), matrix(c(1,rMZ,rMZ,1),2,2))
	DZ <- mvrnorm (1000, c(0,0), matrix(c(1,rDZ,rDZ,1),2,2))

	selVars <- c('t1','t2')
	dimnames(DataMZ) <- list(NULL,selVars)
	dimnames(DataDZ) <- list(NULL,selVars)
	summary(DataMZ)
	summary(DataDZ)
	colMeans(DataMZ,na.rm=TRUE)
	colMeans(DataDZ,na.rm=TRUE)
	cov(DataMZ,use="complete")
	cov(DataDZ,use="complete")

	
We typically start with fitting a saturated model, estimating means, variances and covariances separately by order of the twins (twin 1 vs twin 2) and by zygosity (MZ vs DZ pairs), to establish the likelihood of the data.  This is essentially similar to the optimization script discussed above, except that we now have two variables (same variable for twin 1 and twin 2) and two groups (MZ and DZ). Thus, the saturated model will have two matrices for the expected means of MZs and DZs, and two for the expected covariances, generated from multiplying a lower triangular matrix with its transpose.  The raw data are read in using the ``mxData`` command, and the corresponding objective function ``mxFIMLObjective`` applied.  

.. code-block:: r

	mxModel("MZ",
		mxMatrix(
			type="Full", 
			nrow=1, 
			ncol=2, 
			free=TRUE, 
			values=c(0,0), 
			dimnames=list(NULL,selVars), 
			name="expMeanMZ"), 
		mxMatrix(
			type="Lower", 
			nrow=2, 
			ncol=2,
			free=TRUE
			values=.5,
			dimnames=list(NULL, selVars), 
			name="CholMZ"), 
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

Note that the ``mxModel`` statement for the DZ twins is almost identical to that for MZ twins, except for the names of the objects and data.  If the arguments to the OpenMx command are given in the default order (see i.e. ``?mxMatrix`` to go to the help/reference page for that command), then it is not necessary to include the name of the argument.  Given we skip a few optional arguments, the argument names ``dimnames=`` and ``name=`` are included to refer to the right arguments.  For didactic purposes, we prefer the formatting used for the MZ group, with soft tabs and each argument on a separate line, etc.  (see list of formatting rules).  However, the experienced user may want to use a more compact form, as the one used for the DZ group.

.. code-block:: r            

    mxModel("DZ",
        mxMatrix("Full", 1, 2, T, c(0,0), dimnames=list(NULL, selVars), name="expMeanDZ"), 
        mxMatrix("Lower", 2, 2, T, .5, dimnames=list(NULL, selVars), name="CholDZ"), 
        mxAlgebra(CholDZ %*% t(CholDZ), name="expCovDZ", dimnames=list(selVars, selVars)), 
        mxData(DataDZ, type="raw"), 
        mxFIMLObjective("expCovDZ", "expMeanDZ")),

The two models are then combined in a 'super'model which includes them as arguments.  Additional arguments are an ``mxAlgebra`` statement to add the objective funtions/likelihood of the two submodels.  To evaluate them simultaneously, we use the ``mxAlgebraObjective`` with the previous algebra as its argument.  The ``mxRun`` command is used to start optimization.

.. code-block:: r 

	twinSatModel <- mxModel("twinSat",
		mxModel("MZ",
			mxMatrix("Full", 1, 2, T, c(0,0), dimnames=list(NULL, selVars), name="expMeanMZ"), 
			mxMatrix("Lower", 2, 2, T, .5, dimnames=list(NULL, selVars), name="CholMZ"), 
			mxAlgebra(CholMZ %*% t(CholMZ), name="expCovMZ", dimnames=list(selVars, selVars)), 
			mxData(DataMZ, type="raw"), 
			mxFIMLObjective("expCovMZ", "expMeanMZ")),
		mxModel("DZ",
			mxMatrix("Full", 1, 2, T, c(0,0), dimnames=list(NULL, selVars), name="expMeanDZ"), 
			mxMatrix("Lower", 2, 2, T, .5, dimnames=list(NULL, selVars), name="CholDZ"), 
			mxAlgebra(CholDZ %*% t(CholDZ), name="expCovDZ", dimnames=list(selVars, selVars)), 
			mxData(DataDZ, type="raw"), 
			mxFIMLObjective("expCovDZ", "expMeanDZ")),
		mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
		mxAlgebraObjective("twin"))
	twinSatFit <- mxRun(twinSatModel)

It is always helpful/advised to check the model specifications before interpreting the output.  Here we are interested in the values for the expected mean vectors and covariance matrices, and the goodness-of-fit statistics, including the likelihood, degrees of freedom, and any other derived indices.

.. code-block:: r

    ExpMeanMZ <- mxEval(MZ.expMeanMZ, twinSatFit)
    ExpCovMZ <- mxEval(MZ.expCovMZ, twinSatFit)
    ExpMeanDZ <- mxEval(DZ.expMeanDZ, twinSatFit)
    ExpCovDZ <- mxEval(DZ.expCovDZ, twinSatFit)
    LL_Sat <- mxEval(objective, twinSatFit)

Before we move on to fit the ACE model to the same data, we may want to test some of the assumptions of the twin model, i.e. that the means and variances are the same for twin 1 and twin 2, and that they are the same for MZ and DZ twins.  This can be done as an omnibus test, or stepwise.  Let us start by equating the means for both twins, separately in the two groups.  We accomplish this by using the same label (just one label which will be reused by R) for the two free parameters for the means per group.  As the majority of the previous script stays the same, we start by copying the old model into a new one.  We then include the arguments of the model that require a change.

.. code-block:: r 

    twinSatModelSub1 <- mxModel(twinSatModel,
        mxModel("MZ",
            mxMatrix("Full", 1, 2, T, 0, "mMZ", dimnames=list(NULL, selVars), name="expMeanMZ"), 
        mxModel("DZ", 
            mxMatrix("Full", 1, 2, T, 0, "mDZ", dimnames=list(NULL, selVars), name="expMeanDZ"))
    twinSatFitSub1 <- mxModel(twinSatModelSub1)

If we want to test if we can equate both means and variances across twin order and zygosity at once, we will end up with the following specification.  Note that we use the same label across models for elements that need to be equated.

.. code-block:: r 

    twinSatModelSub2 <- mxModel(twinSatModelSub1,
        mxModel("MZ",
            mxMatrix("Full", 1, 2, T, 0, "mean", dimnames=list(NULL, selVars), name="expMeanMZ"), 
            mxMatrix("Lower", 2, 2, T, .5, labels= c("var","MZcov","var"), 
                dimnames=list(NULL, selVars), name="CholMZ"), 
        mxModel("DZ", 
            mxMatrix("Full", 1, 2, T, 0, "mean", dimnames=list(NULL, selVars), name="expMeanDZ"), 
            mxMatrix("Lower", 2, 2, T, .5, labels= c("var","DZcov","var"), 
                dimnames=list(NULL, selVars), name="CholDZ"))
    twinSatFitSub2 <- mxRun(twinSatModelSub2)

We can compare the likelihood of this submodel to that of the fully saturated model or the previous submodel using the results from ``mxEval`` commands with regular R algebra.  A summary of the model parameters, estimates and goodness-of-fit statistics can also be obtained using ``summary(twinSatFit)``.

.. code-block:: r

	LL_Sat <- mxEval(objective, twinSatFit)
	LL_Sub1 <- mxEval(objective, twinSatFitSub1)
	LRT1= LL_Sub1 - LL_Sat
	LL_Sub2 <- mxEval(objective, twinSatFitSub1)
	LRT2= LL_Sub2 - LL_Sat

Now, we are ready to specify the ACE model to test which sources of variance significantly contribute to the phenotype and estimate their best value.  The structure of this script is going to mimic that of the saturated model.  The main difference is that we no longer estimate the variance-covariance matrix directly, but express it as a function of the three sources of variance, **A**, **C** and **E**.  As the same sources are used for the MZ and the DZ group, the matrices which will represent them are part of the 'super'model.  As these sources are variances, which need to be positive, we typically use a Cholesky decomposition of the standard deviations (and effectively estimate **a** rather then **a^2**, see later for more in depth coverage).  Thus, we specify three separate matrices for the three sources of variance using the ``mxMatrix`` command and 'calculate' the variance components with the ``mxAlgebra`` command.  Note that there are a variety of ways to specify this model, we have picked one that corresponds well to previous Mx code, and has some intuitive appeal.

.. code-block:: r

	#Specify ACE Model
	twinACEModel <- mxModel("twinACE", 
		mxMatrix("Full", 1, 2, T, 20, "mean", dimnames=list(NULL, selVars), name="expMean"), 
		# Matrix expMean for expected mean vector for MZ and DZ twins    
		mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="a", name="X"),
		mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="c", name="Y"),
		mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=.6, label="e", name="Z"),
		# Matrices X, Y, and Z to store the a, c, and e path coefficients
		mxMatrix("Full", nrow=1, ncol=1, free=FALSE, values=.5, name="h"),
		mxAlgebra(X * t(X), name="A"),
		mxAlgebra(Y * t(Y), name="C"),
		mxAlgebra(Z * t(Z), name="E"),
		# Matrixes A, C, and E to compute A, C, and E variance components
		mxAlgebra(rbind(cbind(A+C+E   , A+C),
		cbind(A+C     , A+C+E)), dimnames = list(selVars, selVars), name="expCovMZ"),
		# Matrix expCOVMZ for expected covariance matrix for MZ twins
		mxAlgebra(rbind(cbind(A+C+E   , h%x%A+C),
		cbind(h%x%A+C , A+C+E)), dimnames = list(selVars, selVars), name="expCovDZ"),
		# Matrix expCOVMZ for expected covariance matrix for DZ twins
		mxModel("MZ",
			mxData(DataMZ, type="raw"), 
			mxFIMLObjective("twinACE.expCovMZ", "twinACE.expMean")),
		mxModel("DZ", 
			mxData(DataDZ, type="raw"), 
			mxFIMLObjective("twinACE.expCovDZ", "twinACE.expMean")),
		mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
		mxAlgebraObjective("twin"))
	twinACEFit <- mxRun(twinACEModel)

Relevant output can be generate with ``print`` or ``summary`` statements or specific output can be requested using the ``mxEval`` command.  Typically we would compare this model back to the saturated model to interpret its goodness-of-fit.  Parameter estimates are obtained and can easily be standardized.  A typical analysis would likely include the following output.

.. code-block:: r

	LL_ACE <- mxEval(objective, twinACEFit)
	LRT_ACE= LL_ACE - LL_Sat
	
	#Retrieve expected mean vector and expected covariance matrices
		MZc <- mxEval(expCovMZ, twinACEFit)
		DZc <- mxEval(expCovDZ, twinACEFit)
		M   <- mxEval(expMean, twinACEFit)
	#Retrieve the A, C, and E variance components
		A <- mxEval(A, twinACEFit)
		C <- mxEval(C, twinACEFit)
		E <- mxEval(E, twinACEFit)
	#Calculate standardized variance components
		V <- (A+C+E)
		a2 <- A/V
		c2 <- C/V
		e2 <- E/V
	#Build and print reporting table with row and column names
		ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2)) 
		ACEest <- data.frame(ACEest, row.names=c("Variance Components","Standardized VC"))
		names(ACEest)<-c("A", "C", "E")
	 	ACEest; LL_ACE; LRT_ACE

Similarly to fitting submodels from the saturated model, we typically fit submodels of the ACE model to test the significance of the sources of variance.  One example is testing the significance of shared environmental factors by dropping the free parameter for ``c`` (fixing it to zero).  We call up the previous model and include the new specification for the matrix to be changed, and rerun.

.. code-block:: r

	twinAEModel <- mxModel(twinACEModel,
	    mxMatrix("Full", nrow=1, ncol=1, free=F, values=0, label="c", name="Y"))
	twinAEFit <- mxRun(twinAEModel)

We discuss twin analysis examples in more detail in the example code.
