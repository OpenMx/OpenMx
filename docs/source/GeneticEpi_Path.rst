    .. _geneticepidemiology-path-specification:

Genetic Epidemiology, Path Specification
=========================================

Mx was and OpenMx is probably the most popular statistical modeling package in the behavior genetics field, as it was conceived with genetic models in mind, which rely heavily on multiple groups.  We introduce here an OpenMx script for the basic genetic model in genetic epidemiologic research, the ACE model.  This model assumes that the variability in a phenotype, or observed variable,  can be explained by differences in genetic and environmental factors, with **A** representing additive genetic factors, **C** shared/common environmental factors and **E** unique/specific environmental factors (see Neale & Cardon 1992, for a detailed treatment).  To estimate these three sources of variance, data have to be collected on relatives with different levels of genetic and environmental similarity to provide sufficient information to identify the parameters.  One such design is the classical twin study, which compares the similarity of identical (monozygotic, MZ) and fraternal (dizygotic, DZ) twins to infer the role of **A**, **C** and **E**.

The example starts with the ACE model and includes one submodel, the AE model. It is available in the following file:

* http://openmx.psyc.virginia.edu/svn/trunk/demo/UnivariateTwinAnalysis_PathRaw.R

A parallel version of this example, using matrix specification of models rather than paths, can be found here:

* http://openmx.psyc.virginia.edu/svn/trunk/demo/UnivariateTwinAnalysis_MatrixRaw.R


ACE Model: a Twin Analysis
--------------------------

A twin analysis is a typical example of multiple groups, in this case MZ twins and DZ twins, with different expectations for the covariance structure (and possibly means).  We illustrate the model here with the corresponding two path diagrams:

.. image:: graph/TwinACEModelMZ.png
    :height: 2.5in
    
.. image:: graph/TwinACEModelDZ.png
    :height: 2.5in


Data
^^^^

Let us assume you have collected data on a large sample of twin pairs for your phenotype of interest.  For illustration purposes, we use Australian data on body mass index (BMI) which are saved in a text file 'myTwinData.txt'.  We use R to read the data into a data.frame and define the objects ``selVars`` for the variables selected for analysis, and ``aceVars`` for the latent variables to simplify the OpenMx code.  We then create two subsets of the data for MZ females (mzData) and DZ females (dzData) respectively with the code below, and generate some descriptive statistics, namely means and covariances.

.. code-block:: r

    # Load Data
    data(twinData)

    # Select Variables for Analysis
    selVars   <- c('bmi1','bmi2')
    aceVars   <- c("A1","C1","E1","A2","C2","E2")

    # Select Data for Analysis
    mzData    <- subset(twinData, zyg==1, selVars)
    dzData    <- subset(twinData, zyg==3, selVars)

    # Generate Descriptive Statistics
    colMeans(mzData,na.rm=TRUE)
    colMeans(dzData,na.rm=TRUE)
    cov(mzData,use="complete")
    cov(dzData,use="complete")


Model Specification
^^^^^^^^^^^^^^^^^^^

There are different ways to draw a path diagram of the ACE model.  The most commonly used approach is with the three latent variables in circles at the top, separately for twin 1 and twin 2 respectively called **A1**, **C1**, **E1** and **A2**, **C2**, **E2**.  The latent variables are connected to the observed variables in boxes at the bottom, representing the measures for twin 1 and twin 2: **T1** and **T2**, by single-headed arrows from the latent to the manifest variables.  Path coefficients **a**, **c** and **e** are estimated but constrained to be the same for twin 1 and twin 2, as well as for MZ and DZ twins.  As MZ twins share all their genotypes, the double-headed path connecting **A1** and **A2** is fixed to one in the MZ diagram.  DZ twins share on average half their genes, and as a result the corresponding path is fixed to 0.5 in the DZ diagram.  Environmental factors that are shared between twins are assumed to increase similarity between twins to the same extent in MZ and DZ twins (equal environments assumption), thus the double-headed path connecting **C1** and **C2** is fixed to one in both diagrams above.  The unique environmental factors are by definition uncorrelated between twins.

Let's go through the paths specification step by step.  First, we start with the ``require(OpenMx)`` statement.  We include the full code here.  As MZ and DZ have to be evaluated together, the models for each will be arguments of a bigger model.  Given the diagrams for the MZ and the DZ group look rather similar, we start by specifying all the common elements  which will be stored in a list *paths* and then shared with the two submodels for each of the twin types, defined in separate ``mxModel`` commands.  The latter two ``MxModel`` objects (*modelMZ* and *modelDZ*) are arguments of the overall model, and will be saved together in the R object *modelACE* and thus be run together.

.. code-block:: r

    require(OpenMx)    
    
    # Path objects for Multiple Groups
    manifestVars=selVars
    latentVars=aceVars
    latVariances <- mxPath( from=aceVars, arrows=2, free=FALSE, values=1 )                                             # variances of latent variables
    latMeans     <- mxPath( from="one", to=aceVars, arrows=1, free=FALSE, values=0 )                                   # means of latent variables
    obsMeans     <- mxPath( from="one", to=selVars, arrows=1, free=TRUE, values=20, labels="mean" )                    # means of observed variables
    pathAceT1    <- mxPath( from=c("A1","C1","E1"), to="bmi1", arrows=1, free=TRUE, values=.5,  label=c("a","c","e") ) # path coefficients for twin 1
    pathAceT2    <- mxPath( from=c("A2","C2","E2"), to="bmi2", arrows=1, free=TRUE, values=.5,  label=c("a","c","e") ) # path coefficients for twin 2
    covC1C2      <- mxPath( from="C1", to="C2", arrows=2, free=FALSE, values=1 )                                       # covariance between C1 & C2
    covA1A2_MZ   <- mxPath( from="A1", to="A2", arrows=2, free=FALSE, values=1 )                                       # covariance between A1 & A2 in MZ twins
    covA1A2_DZ   <- mxPath( from="A1", to="A2", arrows=2, free=FALSE, values=.5 )                                      # covariance between A1 & A2 in DZ twins

    # Data objects for Multiple Groups
    dataMZ       <- mxData( observed=mzData, type="raw" )
    dataDZ       <- mxData( observed=dzData, type="raw" )

    # Combine Groups
    paths        <- list( latVariances, latMeans, obsMeans, pathAceT1, pathAceT2, covC1C2 )
    modelMZ      <- mxModel(model="MZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_MZ, dataMZ )
    modelDZ      <- mxModel(model="DZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_DZ, dataDZ )
    minus2ll     <- mxAlgebra( expression=MZ.fitfunction + DZ.fitfunction, name="minus2loglikelihood" )
    obj          <- mxFitFunctionAlgebra( "minus2loglikelihood" )
    modelACE     <- mxModel(model="ACE", modelMZ, modelDZ, minus2ll, obj )

    # Run Model
    fitACE       <- mxRun(modelACE)
    sumACE       <- summary(fitACE)
    

Now we will discuss the script line by line.  For further details on RAM, see [RAM1990].  Each line can be pasted into R, and then evaluated together once the whole model is specified.  Models specifying paths are translated into 'RAM' specifications for optimization, indicated by using the ``type="RAM"`` within the ``mxModel`` statements.  We start the path diagram specification by providing the names for the manifest variables in ``manifestVars`` and the latent variables in ``latentVars``.  We use here the ``selVars`` and ``aceVars`` objects that we created previously when preparing the data.

    ..[RAM1990]  McArdle, J.J. & Boker, S.M. (1990). RAMpath: Path diagram software. Denver: Data Transforms Inc.
    

.. code-block:: r

	        manifestVars=selVars
	        latentVars=aceVars

We start by specifying paths for the variances and means of the latent variables.  These include double-headed arrows from each latent variable back to itself, fixed at one.

.. code-block:: r        

    # variances of latent variables
    latVariances <- mxPath( from=aceVars, arrows=2, free=FALSE, values=1 )

and single-headed arrows from the triangle (with a fixed value of one) to each of the latent variables, fixed at zero. 

.. code-block:: r        

    # means of latent variables
    latMeans     <- mxPath( from="one", to=aceVars, arrows=1, free=FALSE, values=0 )

Next we specify paths for the means of the observed variables using single-headed arrows from ``one`` to each of the manifest variables.  These are set to be free and given a start value of 20.  As we use the same label ("mean") for the two means, they are constrained to be equal.  Remember that R 'recycles'.

.. code-block:: r        

    # means of observed variables
    obsMeans     <- mxPath( from="one", to=selVars, arrows=1, free=TRUE, values=20, labels="mean" )

The main paths of interest are those from each of the latent variables to the respective observed variable.  These are also estimated (thus all are set free), get a start value of 0.5 and appropriate labels.  We chose the start value of .5 by dividing the observed variance, here about .7-.8 in three for the three sources of variance, and then taking the square root as we're estimating the path coefficients, but these are squared to obtain their contribution to the variance.

.. code-block:: r        

    # path coefficients for twin 1
    pathAceT1    <- mxPath( from=c("A1","C1","E1"), to="bmi1", arrows=1, free=TRUE, values=.5,  label=c("a","c","e") )
    # path coefficients for twin 2
    pathAceT2    <- mxPath( from=c("A2","C2","E2"), to="bmi2", arrows=1, free=TRUE, values=.5,  label=c("a","c","e") )
    
As the common environmental factors are by definition the same for both twins, we fix the correlation between **C1** and **C2** to one.    

.. code-block:: r        

    # covariance between C1 & C2
    covC1C2      <- mxPath( from="C1", to="C2", arrows=2, free=FALSE, values=1 )

Next we create the paths that are specific to the MZ group or the DZ group and are later included into the respective models, ``modelMZ`` and ``modelDZ``, which are combined in ``modelACE``.   In the MZ model we add the path for the correlation between **A1** and **A2** which is fixed to one.  In the DZ model the correlation between **A1** and **A2** is fixed to 0.5 instead.

.. code-block:: r

    # covariance between A1 & A2 in MZ's
    covA1A2_MZ   <- mxPath( from="A1", to="A2", arrows=2, free=FALSE, values=1 )
    # covariance between A1 & A2 in DZ's
    covA1A2_DZ   <- mxPath( from="A1", to="A2", arrows=2, free=FALSE, values=.5 )

That concludes the specification of the paths from which the models will be generated for MZ and DZ twins separately.  Next we move to the ``mxData`` commands that call up the data.frame with the MZ raw data, *mzData*, and the DZ raw data, *dzData*, respectively, with the type specified explicitly as *raw*.  These are stored in two MxData objects.

.. code-block:: r

    dataMZ       <- mxData( observed=mzData, type="raw" )
    dataDZ       <- mxData( observed=dzData, type="raw" )

As we indicated earlier, we're collecting all the mxPaths objects that are in common between the two models in a list called *paths*, which will then be included in the respective models that we'll build next with the ``mxModel`` statements.  First we give the model a name, "MZ", to refer back to it later when we need to add the fit functions.  Next we tell OpenMx that we're specifying a path model by using the RAM ``type``, which requires us to include both the ``manifestVars`` and the ``latentVars`` arguments.  Then we include the list of paths generated before that are common between the two models, and the path that is specific to either the MZ or the DZ model.  Last we add the data objects for the MZ and DZ group respectively.

.. code-block:: r    
    
    # Combine Groups
    paths        <- list( latVariances, latMeans, obsMeans, pathAceT1, pathAceT2, covC1C2 )
    modelMZ      <- mxModel(model="MZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_MZ, dataMZ )
    modelDZ      <- mxModel(model="DZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_DZ, dataDZ )

Finally, both models need to be evaluated simultaneously.  We generate the sum of the fit functions for the two groups, using ``mxAlgebra``, and use the result (*minus2loglikelihood*) as argument of the ``mxFitFunctionAlgebra`` command.  We specify a new ``mxModel`` - with a new name using the ``model=""`` notation, which has the *modelMZ* and *modelDZ* as its arguments.  We also include the objects summing the likelihood and evaluating it.

.. code-block:: r        

    minus2ll     <- mxAlgebra( expression=MZ.fitfunction + DZ.fitfunction, name="minus2loglikelihood" )
    obj          <- mxFitFunctionAlgebra( "minus2loglikelihood" )
    modelACE     <- mxModel(model="ACE", modelMZ, modelDZ, minus2ll, obj ) 
    

Model Fitting
^^^^^^^^^^^^^
        
We need to invoke the ``mxRun`` command to start the model evaluation and optimization.  Detailed output will be available in the resulting object, which can be obtained by a ``print()`` statement, or a more succinct output can be obtained with the ``summary`` function.

.. code-block:: r        

    #Run ACE model
    fitACE       <- mxRun(modelACE)
    sumACE       <- summary(fitACE)

Often, however, one is interested in specific parts of the output.  In the case of twin modeling, we typically will inspect the likelihood, the expected covariance matrices and mean vectors, the parameter estimates, and possibly some derived quantities, such as the standardized variance components, obtained by dividing each of the components by the total variance.  Note in the code below that the ``mxEval`` command allows easy extraction of the values in the various matrices which form the first argument, with the model name as second argument.  Once these values have been put in new objects, we can use any regular R expression to derive further quantities or organize them in a convenient format for including in tables.  Note that helper functions could easily (and will likely) be written for standard models to produce 'standard' output. 

.. code-block:: r

    # Generate & Print Output
    A  <- mxEval(a*a, fitACE)                         # additive genetic variance, a^2
    C  <- mxEval(c*c, fitACE)                         # shared environmental variance, c^2
    E  <- mxEval(e*e, fitACE)                         # unique environmental variance, e^2
    V  <- (A+C+E)                                     # total variance
    a2 <- A/V                                         # standardized A
    c2 <- C/V                                         # standardized C
    e2 <- E/V                                         # standardized E
    estACE <- rbind(cbind(A,C,E),cbind(a2,c2,e2))     # table of estimates
    LL_ACE <- mxEval(fitfunction, fitACE)             # likelihood of ACE model

Alternative Models: an AE Model
-------------------------------

To evaluate the significance of each of the model parameters, nested submodels are fit in which the parameters of interest are fixed to zero.  If the likelihood ratio test between the two models (one including the parameter and the other not) is significant, the parameter that is dropped from the model significantly contributes to the variance of the phenotype in question.  Here we show how we can fit the AE model as a submodel with a change in the two ``mxPath`` commands.  We re-specify the path from **C1** to **bmi1** to be fixed to zero, and do the same for the path from **C2** to **bmi2**.  We need to rebuild both modelMZ and  modelDZ, so that they are now built with the changed paths, as well as the overall model which we now call modelAE.  We can run this model in the same way as before, by combining the fit functions of the two groups and generate similar summaries of the results.

.. code-block:: r

    #Run AE model
    pathAceT1    <- mxPath( from=c("A1","C1","E1"), to="bmi1", arrows=1, free=c(T,F,T), values=c(.6,0,.6),  label=c("a","c","e") ) # path coefficients for twin 1
    pathAceT2    <- mxPath( from=c("A2","C2","E2"), to="bmi2", arrows=1, free=c(T,F,T), values=c(.6,0,.6),  label=c("a","c","e") ) # path coefficients for twin 2

    # Combine Groups
    paths        <- list( latVariances, latMeans, obsMeans, pathAceT1, pathAceT2, covC1C2 )
    modelMZ      <- mxModel(model="MZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_MZ, dataMZ )
    modelDZ      <- mxModel(model="DZ", type="RAM", manifestVars=selVars, latentVars=aceVars, paths, covA1A2_DZ, dataDZ )
    modelAE      <- mxModel(model="AE", modelMZ, modelDZ, minus2ll, obj )

    # Run Model
    fitAE        <- mxRun(modelAE)
    sumAE        <- summary(fitAE)

    # Generate & Print Output
    A  <- mxEval(a*a, fitAE)
    C  <- mxEval(c*c, fitAE)
    E  <- mxEval(e*e, fitAE)
    V  <- (A+C+E)
    a2 <- A/V
    c2 <- C/V
    e2 <- E/V
    estAE <- rbind(cbind(A, C, E),cbind(a2, c2, e2))
    LL_AE <- mxEval(fitfunction, fitAE)
    LRT_ACE_AE <- LL_AE - LL_ACE
    estACE
    estAE
    LRT_ACE_AE
    
We use a likelihood ratio test (or take the difference between -2 times the log-likelihoods of the two models, for the difference in degrees of freedom) to determine the best fitting model.  In this example, the Chi-square likelihood ratio test is 0 for 1 degree of freedom, indicating the the *c* parameter does not contribute to the variance at all.  This can also be seen in the 0 estimates for the *c* parameter in the ACE model and identical parameters for *a* and *e* in the ACE and AE models.

While the approach outlined above works just fine, the same can be accomplished with the ``omxSetParameters`` helper function, that allows the user to specify a parameter label in a model whose attributes are changed, in this case by setting ``free`` to FALSE and ``values`` to 0.  Prior to making this changed, we copied the original model into a new model and gave it a new name, so that we have separate model objects for the two nested models that can then be compared with ``mxCompare``.

.. code-block:: r
    
    modelAE    <- mxModel( fitACE, name="AE" )
    modelAE    <- omxSetParameters( modelAE, labels="c", free=FALSE, values=0 )
    fitAE      <- mxRun(modelAE)
    sumAE      <- summary(fitAE)
    mxCompare(fitACE, fitAE)


See :ref:`geneticepidemiology-matrix-specification` for matrix specification of these models.
