Genetic Epidemiology, Matrix Specification
==========================================

Mx is probably most popular in the behavior genetics field, as it was conceived with genetic models in mind, which rely heavily on multiple groups.  We introduce here an OpenMx script for the basic genetic model in genetic epidemiologic research, the ACE model.  This model assumes that the variability in a phenotype, or observed variable, of interest can be explained by differences in genetic and environmental factors, with A representing additive genetic factors, C shared/common environmental factors and E unique/specific environmental factors (see Neale & Cardon 1992, for a detailed treatment).  To estimate these three sources of variance, data have to be collected on relatives with different levels of genetic and environmental similarity to provide sufficient information to identify the parameters.  One such design is the classical twin study, which compares the similarity of identical (monozygotic, MZ) and fraternal (dizygotic, DZ) twins to infer the role of **A**, **C** and **E**.

The example starts with the ACE model and includes one submodel, the AE model. It is available in the following file:

* UnivariateTwinAnalysis_MatrixRaw.R

A parallel version of this example, using path specification of models rather than matrices, can be found here link.

ACE Model: a Twin Analysis
--------------------------

Data
^^^^

Let us assume you have collected data on a large sample of twin pairs for your phenotype of interest.  For illustration purposes, we use Australian data on body mass index (BMI) which are saved in a text file 'myTwinData.txt'.  We use R to read the data into a data.frame and to create two subsets of the data for MZ females (mzfData) and DZ females (dzfData) respectively with the code below.

.. code-block:: r

    require(OpenMx)

    #Prepare Data
    twinData <- read.table("myTwinData.txt", header=T, na.strings=".")
    twinVars <- c('fam','age','zyg','part','wt1','wt2','ht1','ht2','htwt1','htwt2','bmi1','bmi2')
    summary(twinData)
    selVars <- c('bmi1','bmi2')
    mzfData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2)))
    dzfData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2)))


Model Specification
^^^^^^^^^^^^^^^^^^^

There are a variety of ways to set up the ACE model.  The most commonly used approach in Mx is to specify three matrices for each of the three sources of variance.  The matrix **X** represents the additive genetic path 'a', the **Y** matrix is used for the shared environmental path 'c' and the matrix **Z** for the unique environmental path 'e'.  The expected variances and covariances between member of twin pairs are typically expressed in variance components (or the square of the path coefficients, i.e. 'a^2', 'c^2' and 'e^2').  These quantities can be calculated using matrix algebra, by multiplying the **X** matrix by its transpose **t(A)**.  Note that the transpose is not strictly needed in the univariate case, but will allow easier transition to the multivariate case.  We then use matrix algebra again to add the relevant matrices corresponding to the expectations for each of the statistics of the observed covariance matrix.  The R functions 'cbind' and 'rbind' are used to concatenate the resulting matrices in the appropriate way.  Let's go through each of the matrices step by step.  They will all form arguments of the ``mxModel``, specified as follows.  Note that we left the comma's at the end of the lines which are necessary when all the arguments are combined prior to running the model.  Each line can be pasted into R, and then evaluated together once the whole model is specified.

.. code-block:: r

    #Fit ACE Model with RawData and Matrix-style Input
    twinACEModel <- mxModel("twinACE",

As the focus is on individual differences, the model for the means is typically simple.  We can estimate each of the means, in each of the two groups (MZ & DZ) as free parameters.  Alternatively, we can establish whether the means can be equated across order and zygosity by fitting submodels to the saturated model.  In this case, we opted to use one 'grand' mean, obtained by assigning the same label to the two elements of the matrix 'expMeanMZ' and the two elements of the matrix 'expMeanDZ', each of which are 'Full' 1x2 matrices with free parameters and start values of 20.  Note again that ``dimnames`` are required for matrices or algebras that generate the expected mean vectors and expected covariance matrices.

.. code-block:: r

    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=T, 
        values=c(20,20), 
        labels= c("mean","mean"), 
        dimnames=list(NULL, selVars), 
        name="expMeanMZ"
    ), 
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2, 
        free=T, 
        values=c(20,20), 
        labels= c("mean","mean"), 
        dimnames=list(NULL, selVars), 
        name="expMeanDZ"
    ), 

Given the current example is univariate (in the sense that we analyze one variable, even though we have measured it in two members of twin pairs), the matrices for the paths 'a', 'c' and 'e', respectively, **X**, **Y** and **Z** are all 'Full' 1x1 matrices assigned the 'free' status and given a .6 starting value.  We also specify the matrix **h** to have a fixed value of 0.5, necessary for the expectation of DZ twins.  

.. code-block:: r

    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=1, 
        free=TRUE, 
        values=.6, 
        label="a", 
        name="X"
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=1, 
        free=TRUE, 
        values=.6, 
        label="c", 
        name="Y"
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=1, 
        free=TRUE, 
        values=.6, 
        label="e", 
        name="Z"
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=1, 
        free=FALSE, 
        values=.5,  
        name="h"
    ),

While the labels in these matrices are given lower case names, similar to the convention that paths have lower case names, the names for the variance component matrices, obtained from multiplying matrices with their transpose have upper case letters 'A', 'C' and 'E' which are distinct  (as R is case-sensitive).

.. code-block:: r

    mxAlgebra(
        expression=X * t(X), 
        name="A"
    ),
    mxAlgebra(
        expression=Y * t(Y), 
        name="C"
    ),
    mxAlgebra(
        expression=Z * t(Z), 
        name="E"
    ), 
        
Previous Mx users will likely be familiar with the look of the expected covariance matrices for MZ and DZ twin pairs.  These 2x2 matrices are built by horizontal and vertical concatenation of the appropriate matrix expressions for the variance, the MZ and the DZ covariance.  In R, concatenation of matrices is accomplished with the 'rbind' and 'cbind' functions.  Thus to represent the matrices in expression ? in R, we use the following code.

.. math::
   :nowrap:

	\begin{eqnarray*}
   covMZ = \left[ \begin{array}{r} a^2+c^2+e^2, a^2+c^2 \\ a^2+c^2, a^2+c^2+e^2 \\ \end{array} \right]
   & covDZ = \left[ \begin{array}{r} a^2+c^2+e^2, .5a^2+c^2 \\ .5a^2+c^2, a^2+c^2+e^2 \\ \end{array} \right]
	\end{eqnarray*}

.. code-block:: r

    mxAlgebra(
        expression=rbind (cbind(A + C + E, A + C),
                          cbind(A + C    , A + C + E)), 
        dimnames = list(selVars, selVars), 
        name="expCovMZ"
    ),
    mxAlgebra(
        expression=rbind (cbind(A + C + E  , h %x% A + C),
                          cbind(h %x% A + C, A + C + E)), 
        dimnames = list(selVars, selVars), 
        name="expCovDZ"
    ),

As the expected covariance matrices are different for the two groups of twins, we specify two ``mxModel`` commands within the 'twinACE' mxModel command.  They are given a name, and arguments for the data and the objective function to be used to optimize the model.  We have set the model up for raw data, and thus will use the ``mxFIMLObjective`` function to evaluate it.  For each model, the ``mxData`` command calls up the appropriate data, and provides a type, here 'raw', and the ``mxFIMLObjective`` command is given the names corresponding to the respective expected covariance matrices and mean vectors, specified above.

.. code-block:: r

    mxModel("MZ",
        mxData(
            observed=mzfData, 
            type="raw"
        ), 
        mxFIMLObjective(
            covariances="twinACE.expCovMZ", 
            means="twinACE.expMeanMZ"
        )
    ),
    mxModel("DZ", 
        mxData(
            observed=dzfData, 
            type="raw"
        ), 
        mxFIMLObjective(
            covariances="twinACE.expCovDZ", 
            means="twinACE.expMeanDZ"
        )
    ),

Finally, both models need to be evaluated simultaneously.  We first generate the sum of the objective functions for the two groups, using ``mxAlgebra``, and then use that as argument of the ``mxAlgebraObjective`` command.

.. code-block:: r

        mxAlgebra(
            expression=MZ.objective + DZ.objective, 
            name="twin"
        ), 
        mxAlgebraObjective("twin")
        )

Model Fitting
^^^^^^^^^^^^^

We need to invoke the ``mxRun`` command to start the model evaluation and optimization.  Detailed output will be available in the resulting object, which can be obtained by a ``print()`` statement.

.. code-block:: r

    #Run ACE model
    twinACEFit <- mxRun(twinACEModel)

Often, however, one is interested in specific parts of the output.  In the case of twin modeling, we typically will inspect the expected covariance matrices and mean vectors, the parameter estimates, and possibly some derived quantities, such as the standardized variance components, obtained by dividing each of the components by the total variance.  Note in the code below that the ``mxEval`` command allows easy extraction of the values in the various matrices/algebras which form the first argument, with the model name as second argument.  Once these values have been put in new objects, we can use and regular R expression to derive further quantities or organize them in a convenient format for including in tables.  Note that helper functions could (and will likely) easily be written for standard models to produce 'standard' output. 

.. code-block:: r

    MZc <- mxEval(expCovMZ, twinACEFit)
    DZc <- mxEval(expCovDZ, twinACEFit)
    M <- mxEval(expMeanMZ, twinACEFit)
    A <- mxEval(A, twinACEFit)
    C <- mxEval(C, twinACEFit)
    E <- mxEval(E, twinACEFit)
    V <- (A+C+E)
    a2 <- A/V
    c2 <- C/V
    e2 <- E/V
    ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
    LL_ACE <- mxEval(objective, twinACEFit)
    
Alternative Models: an AE Model
-------------------------------

To evaluate the significance of each of the model parameters, nested submodels are fit in which these parameters are fixed to zero.  If the likelihood ratio test between the two models is significant, the parameter that is dropped from the model significantly contributes to the phenotype in question.  Here we show how we can fit the AE model as a submodel with a change in one ``mxmMatrix`` command.  First, we call up the previous 'full' model and save it as a new model 'twinAEModel'.  Next we re-specify the matrix **Y** to be fixed to zero.  We can run this model in the same way as before and generate similar summaries of the results.

.. code-block:: r

    #Run AE model
    twinAEModel <- mxModel(twinACEModel, 
        mxMatrix(
            type="Full", 
            nrow=1, 
            ncol=1, 
            free=F, 
            values=0, 
            label="c", 
            name="Y"
        )
        )
    twinAEFit <- mxRun(twinAEModel)

    MZc <- mxEval(expCovMZ, twinAEFit)
    DZc <- mxEval(expCovDZ, twinAEFit)
    A <- mxEval(A, twinAEFit)
    C <- mxEval(C, twinAEFit)
    E <- mxEval(E, twinAEFit)
    V <- (A+C+E)
    a2 <- A/V
    c2 <- C/V
    e2 <- E/V
    AEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
    LL_AE <- mxEval(objective, twinAEFit)

We use a likelihood ratio test (or take the difference between -2 times the log-likelihoods of the two models) to determine the best fitting model, and print relevant output.

.. code-block:: r

    LRT_ACE_AE <- LL_AE-LL_ACE

    #Print relevant output
    ACEest
    AEest
    LRT_ACE_AE
