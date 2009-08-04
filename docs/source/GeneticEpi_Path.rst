Genetic Epidemiology, Path Specification
=========================================

Mx is probably most popular in the behavior genetics field, as it was conceived with genetic models in mind, which rely heavily on multiple groups.  We introduce here an OpenMx script for the basic genetic model in genetic epidemiologic research, the ACE model.  This model assumes that the variability in a phenotype, or observed variable, of interest can be explained by differences in genetic and environmental factors, with A representing additive genetic factors, C shared/common environmental factors and E unique/specific environmental factors (see Neale & Cardon 1992, for a detailed treatment).  To estimate these three sources of variance, data have to be collected on relatives with different levels of genetic and environmental similarity to provide sufficient information to identify the parameters.  One such design is the classical twin study, which compares the similarity of identical (monozygotic, MZ) and fraternal (dizygotic, DZ) twins to infer the role of **A**, **C** and **E**.

The example starts with the ACE model and includes one submodel, the AE model. It is available in the following file:

* UnivariateTwinAnalysis_PathRaw.R

A parallel version of this example, using matrix specification of models rather than paths, can be found here link.

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
    aceVars <- c("A1","C1","E1","A2","C2","E2")
    mzfData <- as.matrix(subset(twinData, zyg==1, c(bmi1,bmi2)))
    dzfData <- as.matrix(subset(twinData, zyg==3, c(bmi1,bmi2)))

Model Specification
^^^^^^^^^^^^^^^^^^^

There are different ways to draw a path diagram of the ACE model.  The most commonly used approach is with the three latent variables in circles at the top, separately for twin 1 and twin 2 respectively called **A1**, **C1**, **E1** and **A2**, **C2**, **E2**.  The latent variables are connected to the observed variables (in boxes) ***bmi1*** and ***bmi2*** at the bottom by single-headed arrows from the latent to the manifest variables.  Path coefficients **a**, **c** and **e** are estimated but constrained to be the same for twin 1 and twin 2, as well as for MZ and DZ twins.  As MZ twins share all their genotypes, the double-headed path connecting **A1** and **A2** is fixed to one.  DZ twins share on average half their genes, as a result the corresponding path is fixed to 0.5 in the DZ diagram.  As environmental factors that are shared between twins are assumed to increase similarity between twin to the same extent in MZ and DZ twins (equal environments assumption), the double-headed path connecting **C1** and **C2** is fixed to one in both diagrams.  The unique environmental factors are by definition uncorrelated between twins.

.. image:: aceMZ.png
    :height: 280
    
.. image:: aceDZ.png
    :height: 280

Let's go through each of the paths specification step by step.  They will all form arguments of the ``mxModel``, specified as follows.  Given the diagrams for the MZ and the DZ group look rather similar, we start by specifying all the common elements which will then be shared with the two submodels for each of the twin types.  Thus we call the first model 'share'.

.. code-block:: r

    #Fit ACE Model with RawData and Path-style Input
    share <- mxModel("share", 
        type="RAM",

Models specifying paths are translated into 'RAM' specifications for optimization, indicated by using the ``type='RAM'``.  For further details on RAM, see ref.  Note that we left the comma's at the end of the lines which are necessary when all the arguments are combined prior to running the model.  Each line can be pasted into R, and then evaluated together once the whole model is specified.  We start the path diagram specification by providing the names for the manifest variables in ``manifestVars`` and the latent varibles in ``latentVars``.  We use here the 'selVars' and 'aceVars' objects that we created before when preparing the data.

.. code-block:: r

        manifestVars=selVars,
        latentVars=aceVars,

We start by specifying paths for the variances and means of the latent variables.  This includes double-headed arrows from each latent variable back to itself, fixed at one, and single-headed arrows from the triangle (with a fixed value of one) to each of the latent variables, fixed at zero.  Next we specify paths for the means of the observed variables using single-headed arrows from 'one' to each of the manifest variables.  These are set to be free and given a start value of 20.  As we use the same label ("mean") for the two means, they are constrained to be equal.  The main paths of interest are those from each of the latent variables to the respective observed variable.  These are also estimated (thus all are set free), get a start value of .6 and appropriate labels.  As the common environmental factors are by definition the same for both twins, we fix the correlation between **C1** and **C2 to one.

.. code-block:: r        
        
        mxPath(
            from=aceVars, 
            arrows=2, 
            free=FALSE, 
            values=1
        ),
        mxPath(
            from="one", 
            to=aceVars, 
            arrows=1, 
            free=FALSE, 
            values=0
        ),
        mxPath(
            from="one", 
            to=selVars, 
            arrows=1, free=TRUE, 
            values=20, 
            labels= c("mean","mean")
        ),
        mxPath(
            from=c("A1","C1","E1"), 
            to="bmi1", 
            arrows=1, 
            free=TRUE, 
            values=.6, 
            label=c("a","c","e")
        ),
        mxPath(
            from=c("A2","C2","E2"), 
            to="bmi2", 
            arrows=1, 
            free=TRUE, 
            values=.6, 
            label=c("a","c","e")
        ),
        mxPath(
            from="C1", to="C2", 
            arrows=2, 
            free=FALSE, 
            values=1
        )
        )

We add the paths that are specific to the MZ group or the DZ group into the respective submodels which will be combined in 'twinACEModel'.  So we have two ``mxModel`` statement within the "twinACE" model statement.  Each of the two models are based on the previously specified "share" model by including it as its first argument.  Then we add the path for the correlation between **A1** and **A2** which is fixed to one for the MZ group.  That concludes the specification of the model for the MZ's, thus we move to the ``mxData`` command that calls up the data.frame with the MZ raw data, with the type specified explicitly.  Given we use the path specification, the objective function uses RAM, thus ``type='RAM'``.  We also give it the model a name to refer back to it later when we need to add the objective functions.  The ``mxModel`` command for the DZ group is very similar, except that the the correlation between **A1** and **A2** is fixed to 0.5 and the DZ data are read in.

.. code-block:: r

    twinACEModel <- mxModel("twinACE", 
        mxModel(share,
            mxPath(
                from="A1", 
                to="A2", 
                arrows=2, 
                free=FALSE, 
                values=1
            ),
            mxData(
                observed=mzfData, 
                type="raw"), 
            type="RAM", 
            name="MZ"
        ),
        mxModel(share, 
            mxPath(
                from="A1", 
                to="A2", 
                arrows=2, 
                free=FALSE, 
                values=.5
            ),
            mxData(
                observed=dzfData, 
                type="raw"
            ), 
            type="RAM", 
            name="DZ"
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

Often, however, one is interested in specific parts of the output.  In the case of twin modeling, we typically will inspect the expected covariance matrices and mean vectors, the parameter estimates, and possibly some derived quantities, such as the standardized variance components, obtained by dividing each of the components by the total variance.  Note in the code below that the ``mxEvaluate`` command allows easy extraction of the values in the various matrices/algebras which form the first argument, with the model name as second argument.  Once these values have been put in new objects, we can use and regular R expression to derive further quantities or organize them in a convenient format for including in tables.  Note that helper functions could (and will likely) easily be written for standard models to produce 'standard' output. 

.. code-block:: r

    MZc <- mxEvaluate(MZ.covariance, twinACEFit)
    DZc <- mxEvaluate(DZ.covariance, twinACEFit)
    M <- mxEvaluate(MZ.means, twinACEFit)
    A <- mxEvaluate(a*a, twinACEFit)
    C <- mxEvaluate(c*c, twinACEFit)
    E <- mxEvaluate(e*e, twinACEFit)
    V <- (A+C+E)
    a2 <- A/V
    c2 <- C/V
    e2 <- E/V
    ACEest <- rbind(cbind(A,C,E),cbind(a2,c2,e2))
    LL_ACE <- mxEvaluate(objective, twinACEFit)

Alternative Models: an AE Model
-------------------------------

To evaluate the significance of each of the model parameters, nested submodels are fit in which these parameters are fixed to zero.  If the likelihood ratio test between the two models is significant, the parameter that is dropped from the model significantly contributes to the phenotype in question.  Here we show how we can fit the AE model as a submodel with a change in two ``mxPath`` commands.  First, we call up the previous 'full' model and save it as a new model 'twinAEModel'.  Next we re-specify the path from **C1** to **bmi1** to be fixed to zero, and do the same for the path from **C2** to **bmi2**.  We can run this model in the same way as before and generate similar summaries of the results.

.. code-block:: r

    #Run AE model
    twinAEModel <- mxModel(twinACEModel, 
        type="RAM",
        manifestVars=selVars,
        latentVars=aceVars,
        mxPath(
            from=c("A1","C1","E1"), 
            to="bmi1", 
            arrows=1, 
            free=c(T,F,T), 
            values=c(.6,0,.6), 
            label=c("a","c","e")
        ),
        mxPath(
            from=c("A2","C2","E2"), 
            to="bmi2", 
            arrows=1, 
            free=c(T,F,T), 
            values=c(.6,0,.6), 
            label=c("a","c","e")
        )
    )
    twinAEFit <- mxRun(twinAEModel)

    MZc <- mxEvaluate(MZ.covariance, twinAEFit)
    DZc <- mxEvaluate(DZ.covariance, twinAEFit)
    M <- mxEvaluate(MZ.means, twinAEFit)
    A <- mxEvaluate(a*a, twinAEFit)
    C <- mxEvaluate(c*c, twinAEFit)
    E <- mxEvaluate(e*e, twinAEFit)
    V <- (A + C + E)
    a2 <- A / V
    c2 <- C / V
    e2 <- E / V
    AEest <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
    LL_AE <- mxEvaluate(objective, twinAEFit)

We use a likelihood ratio test (or take the difference between -2 times the log-likelihoods of the two models) to determine the best fitting model, and print relevant output.

.. code-block:: r

    LRT_ACE_AE <- LL_AE - LL_ACE

    #Print relevant output
    ACEest
    AEest
    LRT_ACE_AE
