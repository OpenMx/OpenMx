.. _factoranalysis-path-specification:

Factor Analysis, Path Specification
=====================================

This example will demonstrate latent variable modeling via the common factor model using path-centric model specification. We'll walk through two applications of this approach: one with a single latent variable, and one with two latent variables. As with previous examples, these two applications are split into four files, with each application represented separately with raw and covariance data. These examples can be found in the following files:

* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/OneFactorModel_PathCov.R
* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/OneFactorModel_PathRaw.R
* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/TwoFactorModel_PathCov.R
* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/TwoFactorModel_PathRaw.R

Parallel versions of this example, using matrix specification of models rather than paths, can be found here:

* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/OneFactorModel_MatrixCov.R
* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/OneFactorModel_MatrixRaw.R
* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/TwoFactorModel_MatrixCov.R
* http://openmx.psyc.virginia.edu/docs/OpenMx/latest/_static/demo/TwoFactorModel_MatrixRaw.R

Common Factor Model
-------------------

The common factor model is a method for modeling the relationships between observed variables believed to measure or indicate the same latent variable. While there are a number of exploratory approaches to extracting latent factor(s), this example uses structural modeling to fit confirmatory factor models. The model for any person and path diagram of the common factor model for a set of variables :math:`x_{1}`-:math:`x_{6}` are given below.

.. math::
   :nowrap:
   
   \begin{eqnarray*} 
   x_{ij} = \mu_{j} + \lambda_{j} * \eta_{i} + \epsilon_{ij}
   \end{eqnarray*}

.. image:: graph/OneFactorModel.png
    :height: 2in

While 19 parameters are displayed in the equation and path diagram above (six manifest variances, six manifest means, six factor loadings and one factor variance), we must constrain either the factor variance or one factor loading to a constant to identify the model and scale the latent variable. As such, this model contains 18 parameters. Unlike the manifest variable examples we've run up until now, this model is not fully saturated. The means and covariance matrix for six observed variables contain 27 degrees of freedom, and thus our model contains 9 degrees of freedom. 

Data
^^^^

Our first step to running this model is to include the data to be analyzed. The data for this example contain nine variables. We'll select the six we want for this model using the selection operators used in previous examples. Both raw and covariance data are included below, but only one is required for any model.

.. cssclass:: input
..
   
.. code-block:: r

    data(myFADataRaw)
    names(myFADataRaw)

    oneFactorRaw <- myFADataRaw[,c("x1", "x2", "x3", "x4", "x5", "x6")]

    myFADataCov <- matrix(
        c(0.997, 0.642, 0.611, 0.672, 0.637, 0.677, 0.342, 0.299, 0.337,
          0.642, 1.025, 0.608, 0.668, 0.643, 0.676, 0.273, 0.282, 0.287,
          0.611, 0.608, 0.984, 0.633, 0.657, 0.626, 0.286, 0.287, 0.264,
          0.672, 0.668, 0.633, 1.003, 0.676, 0.665, 0.330, 0.290, 0.274,
          0.637, 0.643, 0.657, 0.676, 1.028, 0.654, 0.328, 0.317, 0.331,
          0.677, 0.676, 0.626, 0.665, 0.654, 1.020, 0.323, 0.341, 0.349,
          0.342, 0.273, 0.286, 0.330, 0.328, 0.323, 0.993, 0.472, 0.467,
          0.299, 0.282, 0.287, 0.290, 0.317, 0.341, 0.472, 0.978, 0.507,
          0.337, 0.287, 0.264, 0.274, 0.331, 0.349, 0.467, 0.507, 1.059), nrow=9,
    dimnames=list( c("x1","x2","x3","x4","x5","x6","y1","y2","y3"),
                   c("x1","x2","x3","x4","x5","x6","y1","y2","y3")) )

    oneFactorCov <- myFADataCov[c("x1","x2","x3","x4","x5","x6"),
                                c("x1","x2","x3","x4","x5","x6")]

    myFADataMeans <- c(2.988, 3.011, 2.986, 3.053, 3.016, 3.010, 2.955, 2.956, 2.967)
    names(myFADataMeans) <- c("x1","x2","x3","x4","x5","x6","y1","y2","y3")

    oneFactorMeans <- myFADataMeans[1:6]

Model Specification
^^^^^^^^^^^^^^^^^^^

Creating a path-centric factor model will use many of the same functions and arguments used in previous path-centric examples. However, the inclusion of latent variables adds a few extra pieces to our model. Before running a model, the OpenMx library must be loaded into R using either the ``require()`` or ``library()`` function. All objects required for estimation (data, paths, and a model type) are included in their own arguments or functions. This code uses the ``mxModel`` function to create an ``MxModel`` object, which we will then run.

.. cssclass:: input
..
   
.. code-block:: r

    require(OpenMx)

    dataRaw      <- mxData( observed=myFADataRaw, type="raw" )
    # residual variances
    resVars      <- mxPath( from=c("x1","x2","x3","x4","x5","x6"), arrows=2,
                            free=TRUE, values=c(1,1,1,1,1,1),
                            labels=c("e1","e2","e3","e4","e5","e6") ) 
    # latent variance
    latVar       <- mxPath( from="F1", arrows=2,
                            free=TRUE, values=1, labels ="varF1" )
    # factor loadings	
    facLoads     <- mxPath( from="F1", to=c("x1","x2","x3","x4","x5","x6"), arrows=1,
                            free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=c(1,1,1,1,1,1),
                            labels =c("l1","l2","l3","l4","l5","l6") )
    # means
    means        <- mxPath( from="one", to=c("x1","x2","x3","x4","x5","x6","F1"), arrows=1,
                            free=c(T,T,T,T,T,T,FALSE), values=c(1,1,1,1,1,1,0),
                            labels =c("meanx1","meanx2","meanx3",
                                      "meanx4","meanx5","meanx6",NA) ) 

    oneFactorModel <- mxModel("Common Factor Model Path Specification", type="RAM",
                            manifestVars=c("x1","x2","x3","x4","x5","x6"), latentVars="F1",
                            dataRaw, resVars, latVar, facLoads, means)
    

As with previous examples, this model begins with a name ("Common Factor Model Path Specification") for the model and a ``type="RAM"`` argument. The name for the model may be omitted, or may be specified in any other place in the model using the ``name`` argument. Including ``type="RAM"`` allows the ``mxModel`` function to interpret the ``mxPath`` functions that follow and turn those paths into an expected covariance matrix and means vector for the ensuing data. The ``mxData`` function works just as in previous examples, and the following raw data specification is included in the code: 

.. cssclass:: input
..
   
.. code-block:: r

    dataRaw      <- mxData( observed=myFADataRaw, type="raw" )

can be replaced with a covariance matrix and means, like so:

.. cssclass:: input
..
   
.. code-block:: r

    dataCov      <- mxData( observed=oneFactorCov, type="cov", numObs=500,
                            means=oneFactorMeans )
          
The first departure from our previous examples can be found in the addition of the ``latentVars`` argument after the ``manifestVars`` argument. The ``manifestVars`` argument includes the six variables in our observed data. The ``latentVars`` argument provides names for the latent variables (here just one), so that it may be referenced in ``mxPath`` functions.

.. cssclass:: input
..
   
.. code-block:: r

    manifestVars=c("x1","x2","x3","x4","x5","x6")
    latentVars="F1"

Our model is defined by four ``mxPath`` functions. The first defines the residual variance terms for our six observed variables. The ``to`` argument is not required, as we are specifiying two headed arrows both from and to the same variables, as specified in the ``from`` argument. These six variances are all freely estimated, have starting values of 1, and are labeled ``e1`` through ``e6``.

.. cssclass:: input
..
   
.. code-block:: r

    # residual variances
    resVars      <- mxPath( from=c("x1","x2","x3","x4","x5","x6"), arrows=2,
                            free=TRUE, values=c(1,1,1,1,1,1),
                            labels=c("e1","e2","e3","e4","e5","e6") ) 
      
We also must specify the variance of our latent variable. This code is identical to our residual variance code above, with the latent variable ``"F1"`` replacing our six manifest variables.   Alternatively, both could be combined.
      
.. cssclass:: input
..
   
.. code-block:: r

    # latent variance
    latVar       <- mxPath( from="F1", arrows=2,
                            free=TRUE, values=1, labels ="varF1" )
          
Next come the factor loadings. These are specified as asymmetric paths (regressions) of the manifest variables on the latent variable ``"F1"``. As we have to scale the latent variable, the first factor loading has been given a fixed value of one by setting the first elements of the ``free`` and ``values`` arguments to ``FALSE`` and ``1``, respectively. Alternatively, the latent variable could have been scaled by fixing the factor variance to 1 in the previous ``mxPath`` function and freely estimating all factor loadings. The five factor loadings that are freely estimated are all given starting values of 1 and labels ``l2`` through ``l6``.   
          
.. cssclass:: input
..
   
.. code-block:: r

    # factor loadings
    facLoads     <- mxPath( from="F1", to=c("x1","x2","x3","x4","x5","x6"), arrows=1,
                            free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=c(1,1,1,1,1,1),
                            labels =c("l1","l2","l3","l4","l5","l6") )

Lastly, we must specify the mean structure for this model. As there are a total of seven variables in this model (six manifest and one latent), we have the potential for seven means. However, we must constrain at least one mean to a constant value, as there is not sufficient information to yield seven mean and intercept estimates from the six observed means. The six observed variables receive freely estimated intercepts, while the factor mean is fixed to a value of zero in the code below.
     
.. cssclass:: input
..
   
.. code-block:: r

    # means
    means        <- mxPath( from="one", to=c("x1","x2","x3","x4","x5","x6","F1"), arrows=1,
                            free=c(T,T,T,T,T,T,FALSE), values=c(1,1,1,1,1,1,0),
                            labels =c("meanx1","meanx2","meanx3","meanx4","meanx5","meanx6",NA) )

The model can now be run using the ``mxRun`` function, and the output of the model can be accessed from the ``output`` slot of the resulting model.
A summary of the output can be reached using ``summary()``.

.. cssclass:: input
..
   
.. code-block:: r

    oneFactorFit <- mxRun(oneFactorModel)

    oneFactorFit$output
    summary(oneFactorFit)

Two Factor Model
-------------------

The common factor model can be extended to include multiple latent variables. The model for any person and path diagram of the common factor model for a set of variables :math:`x_{1}`-:math:`x_{3}` and :math:`y_{1}`-:math:`y_{3}` are given below.

.. math::
   :nowrap:
   
   \begin{eqnarray*} 
   x_{ij} = \mu_{j} + \lambda_{j} * \eta_{1i} + \epsilon_{ij}\\
   y_{ij} = \mu_{j} + \lambda_{j} * \eta_{2i} + \epsilon_{ij}
   \end{eqnarray*}

.. image:: graph/TwoFactorModel.png
    :height: 2in

Our model contains 21 parameters (six manifest variances, six manifest means, six factor loadings, two factor variances and one factor covariance), but each factor requires one identification constraint. Like in the common factor model above, we will constrain one factor loading for each factor to a value of one. As such, this model contains 19 parameters. The means and covariance matrix for six observed variables contain 27 degrees of freedom, and thus our model contains 8 degrees of freedom. 

The data for the two factor model can be found in the ``myFAData`` files introduced in the common factor model. For this model, we will select three *x* variables (``x1-x3``) and three *y* variables (``y1-y3``).

.. cssclass:: input
..
   
.. code-block:: r

    twoFactorRaw <- myFADataRaw[,c("x1","x2","x3","y1","y2","y3")]

    twoFactorCov <- myFADataCov[c("x1","x2","x3","y1","y2","y3"),
                                c("x1","x2","x3","y1","y2","y3")]

    twoFactorMeans <- myFADataMeans[c(1:3,7:9)]

Specifying the two factor model is virtually identical to the single factor case. The last three variables of our ``manifestVars`` argument have changed from ``"x4","x5","x6"`` to ``"y1","y2","y3"``, which is carried through references to the variables in later ``mxPath`` functions.
 
.. cssclass:: input
..
   
.. code-block:: r

    dataRaw      <- mxData( observed=twoFactorRaw, type="raw" )
    # residual variances
    resVars      <- mxPath( from=c("x1", "x2", "x3", "y1", "y2", "y3"), arrows=2,
                            free=TRUE, values=c(1,1,1,1,1,1),
                            labels=c("e1","e2","e3","e4","e5","e6") ) 
    # latent variances and covariance
    latVars      <- mxPath( from=c("F1","F2"), arrows=2, connect="unique.pairs",
                            free=TRUE, values=c(1,.5,1), labels=c("varF1","cov","varF2") )
    # factor loadings for x variables	
    facLoadsX    <- mxPath( from="F1", to=c("x1","x2","x3"), arrows=1,
                            free=c(F,T,T), values=c(1,1,1), labels=c("l1","l2","l3") )
    # factor loadings for y variables
    facLoadsY    <- mxPath( from="F2", to=c("y1","y2","y3"), arrows=1,
                            free=c(F,T,T), values=c(1,1,1), labels=c("l4","l5","l6") )
    # means
    means        <- mxPath( from="one", to=c("x1","x2","x3","y1","y2","y3","F1","F2"), 
                            arrows=1,
                            free=c(T,T,T,T,T,T,F,F), values=c(1,1,1,1,1,1,0,0),
                            labels=c("meanx1","meanx2","meanx3",
                                     "meany1","meany2","meany3",NA,NA) ) 

    twoFactorModel <- mxModel("Two Factor Model Path Specification", type="RAM",
                            manifestVars=c("x1", "x2", "x3", "y1", "y2", "y3"), 
                            latentVars=c("F1","F2"),
                            dataRaw, resVars, latVars, facLoadsX, facLoadsY, means)
  
We've covered the ``type`` argument, ``mxData`` function and ``manifestVars`` and ``latentVars`` arguments previously, so now we will focus on the changes this model makes to the ``mxPath`` functions. The first and last ``mxPath`` functions, which detail residual variances and intercepts, accomodate the changes in manifest and latent variables but carry out identical functions to the common factor model.

.. cssclass:: input
..
   
.. code-block:: r 

    # residual variances
    resVars      <- mxPath( from=c("x1", "x2", "x3", "y1", "y2", "y3"), arrows=2,
                            free=TRUE, values=c(1,1,1,1,1,1),
                            labels=c("e1","e2","e3","e4","e5","e6") ) 
    # means
    means        <- mxPath( from="one", to=c("x1","x2","x3","y1","y2","y3","F1","F2"), 
                            arrows=1,
                            free=c(T,T,T,T,T,T,F,F), values=c(1,1,1,1,1,1,0,0),
                            labels=c("meanx1","meanx2","meanx3",
                                     "meany1","meany2","meany3",NA,NA) )
  
The second, third and fourth ``mxPath`` functions provide some changes to the model. The second ``mxPath`` function specifies the variances and covariance of the two latent variables. Like previous examples, we've omitted the ``to`` argument for this set of two-headed paths. Unlike previous examples, we've set the ``connect`` argument to ``unique.pairs``, which creates all unique paths between the variables. As omitting the ``to`` argument is identical to putting identical variables in the ``from`` and ``to`` arguments, we are creating all unique paths from and to our two latent variables. This results in three paths: from F1 to F1 (the variance of F1), from F1 to F2 (the covariance of the latent variables), and from F2 to F2 (the variance of F2). 

.. cssclass:: input
..
   
.. code-block:: r 

    # latent variances and covariance
    latVars      <- mxPath( from=c("F1","F2"), arrows=2, connect="unique.pairs",
                            free=TRUE, values=c(1,.5,1), labels=c("varF1","cov","varF2") )

  
The third and fourth ``mxPath`` functions define the factor loadings for each of the latent variables. We've split these loadings into two functions, one for each latent variable. The first loading for each latent variable is fixed to a value of one, just as in the previous example.

.. cssclass:: input
..
   
.. code-block:: r 

    # factor loadings for x variables
    facLoadsX    <- mxPath( from="F1", to=c("x1","x2","x3"), arrows=1,
                            free=c(F,T,T), values=c(1,1,1), labels=c("l1","l2","l3") )
    # factor loadings for y variables
    facLoadsY    <- mxPath( from="F2", to=c("y1","y2","y3"), arrows=1,
                            free=c(F,T,T), values=c(1,1,1), labels=c("l4","l5","l6") )

  
The model can now be run using the ``mxRun`` function, and the output of the model can be accessed from the ``$output`` slot of the resulting model. A summary of the output can be reached using ``summary()``.

.. cssclass:: input
..
   
.. code-block:: r

    oneFactorFit <- mxRun(oneFactorModel)

    oneFactorFit$output
    summary(oneFactorFit)

These models may also be specified using matrices instead of paths. See :ref:`factoranalysis-matrix-specification` for matrix specification of these models.
