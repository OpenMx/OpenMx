Factor Models, Path Specification
=====================================

This example will demonstrate latent variable modeling via the common factor model using path-centric model specification. We'll walk through two applications of this approach: one with a single latent variable, and one with two latent variables. As with previous examples, these two applications are split into four files, with each application represented separately with raw and covariance data. These examples can be found in the following files:

* OneFactorModel_PathCov.R
* OneFactorModel_PathRaw.R
* TwoFactorModel_PathCov.R
* TwoFactorModel_PathRaw.R

A parallel version of this example, using matrix specification of models rather than paths, can be found here link.

Common Factor Model
-------------------

The common factor model is a method for modeling the relationships between observed variables believed to measure or indicate the same latent variable. While there are a number of exploratory approaches to extracting latent factor(s), this example uses structural modeling to fit confirmatory factor models. The model for any person and path diagram of the common factor model for a set of variables :math:`x_{1}'-:math:`x_{6}' are given below.

.. math::
   :nowrap:
   
   \begin{eqnarray*} 
   x_{ij} = \mu_{j} + \lambda_{j} * \eta_{i} + \epsilon_{ij}
   \end{eqnarray*}

.. image:: Factor1.png

While 19 parameters are displayed in the equation and path diagram above (6 manifest variances, six manifest means, six factor loadings and one factor variance), we must constrain either the factor variance or one factor loading to a constant to identify the model and scale the latent variable. As such, this model contains 18 parameters. Unlike the manifest variable examples we've run up until now, this model is not fully saturated. The means and covariance matrix for six observed variables contain 27 degrees of freedom, and thus our model contains 9 degrees of freedom. 

Data
----

Our first step to running this model is to put include the data to be analyzed. The data for this example contain nine variables. We'll select the six we want for this model using the selection operators used in previous examples. Both raw and covariance data are included below, but only one is required for any model.

.. code-block:: r

  myFADataRaw <- read.table("myFAData.txt",header=TRUE)

  > names(myFADataRaw)
  [1] "x1" "x2" "x3" "x4" "x5" "x6" "y1" "y2" "y3"

  oneFactorRaw <- myFADataRaw[,c("x1", "x2", "x3", "x4", "x5", "x6"]

  myFADataCov <- matrix(
      c(0.997, 0.642, 0.611, 0.672, 0.637, 0.677, 0.342, 0.299, 0.337,
        0.642, 1.025, 0.608, 0.668, 0.643, 0.676, 0.273, 0.282, 0.287,
        0.611, 0.608, 0.984, 0.633, 0.657, 0.626, 0.286, 0.287, 0.264,
        0.672, 0.668, 0.633, 1.003, 0.676, 0.665, 0.330, 0.290, 0.274,
        0.637, 0.643, 0.657, 0.676, 1.028, 0.654, 0.328, 0.317, 0.331,
        0.677, 0.676, 0.626, 0.665, 0.654, 1.020, 0.323, 0.341, 0.349,
        0.342, 0.273, 0.286, 0.330, 0.328, 0.323, 0.993, 0.472, 0.467,
        0.299, 0.282, 0.287, 0.290, 0.317, 0.341, 0.472, 0.978, 0.507,
        0.337, 0.287, 0.264, 0.274, 0.331, 0.349, 0.467, 0.507, 1.059),
      nrow=9,
      dimnames=list(
          c("x1", "x2", "x3", "x4", "x5", "x6", "y1", "y2", "y3"),
          c("x1", "x2", "x3", "x4", "x5", "x6", "y1", "y2", "y3")),
      )

  oneFactorCov <- myFADataCov[c("x1", "x2", "x3", "x4", "x5", "x6"),c("x1", "x2", "x3", "x4", "x5", "x6")]
  
  myFADataMeans <- c(2.988, 3.011, 2.986, 3.053, 3.016, 3.010, 2.955, 2.956, 2.967)
  
  oneFactorMeans <- myFADataMeans[1:6]

Specifying the Model
--------------------

Creating a path-centric factor model will use many of the same functions and arguments used in previous path-centric examples. However, the inclusion of latent variables adds a few extra pieces to our model. Before running a model, the OpenMx library must be loaded into R using either the ``require()`` or ``library()`` function. All objects required for estimation (data, paths, and a model type) are included in their own arguments or functions. This code uses the ``mxModel`` function to create an ``MxModel`` object, which we'll then run.

.. code-block:: r

  oneFactorModel<-mxModel("Common Factor Model - Path", 
      type="RAM",
      mxData(
          data=oneFactorRaw,
          type="raw"),
      manifestVars=c("x1","x2","x3","x4","x5","x6"),
      latentVars="F1",
      # residual variances
      mxPath(from=c("x1","x2","x3","x4","x5","x6"),
          arrows=2,
          free=TRUE,
          values=c(1,1,1,1,1,1),
          labels=c("e1","e2","e3","e4","e5","e6")
          ),
      # latent variance
      mxPath(from="F1",
          arrows=2,
          free=TRUE,
          values=1,
          labels ="varF1"
          ),
      # factor loadings
      mxPath(from="F1",
          to=c("x1","x2","x3","x4","x5","x6"),
          arrows=1,
          free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
          values=c(1,1,1,1,1,1),
          labels =c("l1","l2","l3","l4","l5","l6")
          ),
      # means
      mxPath(from="one",
          to=c("x1","x2","x3","x4","x5","x6","F1"),
          arrows=1,
          free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE),
          values=c(1,1,1,1,1,1,0),
          labels =c("meanx1","meanx2","meanx3",
              "meanx4","meanx5","meanx6",
              NA)
          )
      ) # close model

As with previous examples, this model begins with a name for the model and a ``type="RAM"`` argument. The name for the model may be omitted, or may be specified an any other place in the model using the ``name`` argument. Including ``type="RAM"`` allows the ``mxModel`` function to interpret the ``mxPath`` functions that follow and turn those paths into an expected covariance matrix and means vector for the ensuing data. The ``mxData`` function works just as in previous examples, and the raw data specification included in the code: 

.. code-block:: r

      mxData(
          data=oneFactorRaw,
          type="raw")
          
can be replaced with a covariance matrix and means, like so:

.. code-block:: r

  oneFactorModel<-mxModel("Common Factor Model - Path", 
      type="RAM",
      mxData(
          data=oneFactorCov,
          type="cov",
          numObs=500,
          means=oneFactorMeans)
          
The first departure from our previous examples can be found in the addition of the ``latentVars`` argument after the ``manifestVars`` argument.









Two Factor Model
-------------------

The common factor model can be extended to include multiple latent variables. The model for any person and path diagram of the common factor model for a set of variables :math:`x_{1}'-:math:`x_{3}' and :math:`y_{1}'-:math:`y_{3}' are given below.

.. math::
   :nowrap:
   
   \begin{eqnarray*} 
   x_{ij} = \mu_{j} + \lambda_{j} * \eta_{1i} + \epsilon_{ij}\\
   y_{ij} = \mu_{j} + \lambda_{j} * \eta_{2i} + \epsilon_{ij}
   \end{eqnarray*}

.. image:: Factor1.png