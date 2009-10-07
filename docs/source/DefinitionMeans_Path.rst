Definition Variables, Path Specification
========================================

This example will demonstrate the use of OpenMx definition variables with the analysis of a simple two group dataset.  What are definition variables?  Essentially, definition variables can be thought of as observed variables that are used to change the statistical model on an individual case basis.  In essence, it is as though one or more variables in the raw data vectors are used to specify the statistical model for that individual.  Many different types of statistical model can be specified in this fashion; some  are readily specified in standard fashion, and some cannot.  To illustrate, we implement a two-group model.  The groups differ in their means but not in their variances and covariances.  This situation could easily be modeled in a regular multiple group fashion - it is only implemented using definition variables to illustrate their use.  The results are verified using summary statistics and an Mx 1.0 script for comparison is also available.

Mean Differences
----------------

The scripts are presented here

* DefinitionMeans_PathRaw.R
* DefinitionMeans_PathRaw.mx

Statistical Model
^^^^^^^^^^^^^^^^^

Algebraically, we are going to fit the following model to the observed x and y variables:

.. math::
   :nowrap:
   
   \begin{eqnarray*} 
   x_{i} = \mu_{x} + \beta_x * def + \epsilon_{xi}\\
   y_{i} = \mu_{y} + \beta_y * def + \epsilon_{yi}
   \end{eqnarray*}

where the residual sources of variance, :math:`\epsilon_{xi}` and :math:`\epsilon_{yi}` covary to the extent :math:`\rho`.  So, the task is to estimate: the two means :math:`\mu_{x}` and :math:`\mu_{y}`; the deviations from these means due to belonging to the group identified by having def set to 1 (as opposed to zero), :math:`\beta_{x}` and :math:`\beta_{y}`; and the parameters of the variance covariance matrix: cov(:math:`\epsilon_{x},\epsilon_{y}`).

Our task is to implement the model shown in the Figure below:

.. image:: def.png
    :height: 400

Data Simulation
^^^^^^^^^^^^^^^

Our first step to running this model is to simulate the data to be analyzed. Each individual is measured on two observed variables, x and y, and a third variable "def" which denotes their group membership with a 1 or a 0.  These values for group membership are not accidental, and must be adhered to in order to obtain readily interpretable results.  Other values such as 1 and 2 would yield the same model fit, but would make the interpretation more difficult.  

.. code-block:: r

    library(MASS)  # to get hold of mvrnorm function 

    set.seed(200)  # to make the simulation repeatable
    n = 500        # sample size, per group
  
    Sigma <- matrix(c(1,.5,.5,1),2,2)
    group1<-mvrnorm(n=n, c(1,2), Sigma)
    group2<-mvrnorm(n=n, c(0,0), Sigma)

We make use of the superb R function ``mvrnorm`` in order to simulate n=500 records of data for each group.  These observations correlate .5 and have a variance of 1, per the matrix Sigma.  The means of x and y in group 1 are 1.0 and 2.0, respectively; those in group 2 are both zero.  The output of the ``mvrnorm`` function calls are matrices with 500 rows and 3 columns, which are stored in group 1 and group 2.  Now we create the definition variable

.. code-block:: r

    # Put the two groups together, create a definition variable, 
    # and make a list of which variables are to be analyzed (selvars)
    y<-rbind(group1,group2)
    dimnames(y)[2]<-list(c("x","y"))
    def<-rep(c(1,0),each=n)
    selvars<-c("x","y")

The objects y and def might be combined in a data frame.  However, in this case we won't bother to do it externally, and simply paste them together in the mxData function call.

Model Specification
^^^^^^^^^^^^^^^^^^^


Before specifying a model, the OpenMx library must be loaded into R using either the ``require()`` or ``library()`` function. This code uses the ``mxModel`` function to create an ``mxModel`` object, which we'll then run.  Note that all the objects required for estimation (data, matrices, and an objective function) are declared within the ``mxModel`` function.  This type of code structure is recommended for OpenMx scripts generally.

.. code-block:: r

    require(OpenMx)
    defmeansmodel<-mxModel("Definition Means via Paths", 
        type="RAM",

The first argument in an ``mxModel`` function has a special function. If an object or variable containing an ``MxModel`` object is placed here, then ``mxModel`` adds to or removes pieces from that model. If a character string (as indicated by double quotes) is placed first, then that becomes the name of the model. Models may also be named by including a ``name`` argument. This model is named ``"DefinitionMeans"``.

The second line of the mxModel function call declares that we are going to be using RAM specification of the model, using directional and bidirectional
path coefficients between the variables. Next, we declare where the data are, and their type, by creating an ``MxData`` object with the ``mxData``
function. This code first references the object where our data are, then uses the ``type`` argument to specify that this is raw data. Analyses using
definition variables have to use raw data, so that the model can be specified on an individual data vector level.

.. code-block:: r

    mxData(
        observed=data.frame(y,def), 
        type="raw"),
    manifestVars=c("x","y"),
    latentVars="DefDummy",

Model specification is carried out using two lists of variables, ``manifestVars`` and ``latentVars``.  Then ``mxPath`` functions are used to specify paths between them. In the present case, we need four mxPath commands to specify the model.  The first is for the variances of the x and y variables, and the second specifies their covariance.  The third specifies a path from the mean vector, always known by the special keword "one", to each of the observed variables, and to the single latent variable "DefDummy".  This last path is specified to contain the definition variable, by virtue of the "data.def" label.  Finally, two paths are specified from the "DefDummy" latent variable to the observed variables.  These parameters estimate the deviation of the mean of those with a data.def value of 1 from that of those with data.def values of zero.

.. code-block:: r

    mxPath(from=c("x","y"), 
        arrows=2,
        free=TRUE,
        values=c(1,.1,1),
        labels=c("Varx","Vary")
    ), # variances
    mxPath(from="x", to="y",
        arrows=2,
        free=TRUE,
        values=c(.1),
        labels=c("Covxy")
    ), # covariances
    mxPath(from="one",
        to=c("x","y","DefDummy"),
        arrows=1,
        free=c(TRUE,TRUE,FALSE),
        values=c(1,1,1),
        labels =c("meanx","meany","data.def")
    ), # means
    mxPath(from="DefDummy",
        to=c("x","y"),
        arrows=1,
        free=c(TRUE,TRUE),
        values=c(1,1),
        labels =c("beta_1","beta_2")
    )) # moderator paths

We can then run the model and examine the output with a few simple commands.

Model Fitting
^^^^^^^^^^^^^^

.. code-block:: r

    # Run the model
    defMeansFit<-mxRun(defMeansModel)
    defMeansFit@matrices

The R object ``defmeansresult`` contains matrices and algebras; here we are interested in the matrices, which can be seen with the ``defmeansresult@matrices`` entry.  In path notation, the unidirectional, one-headed arrows appear in the matrix A, the two-headed arrows in S, and the mean vector single headed arrows in M.

.. code-block:: r

    # Compare OpenMx estimates to summary statistics from raw data, 
    # remembering to knock off 1 and 2 from group 1's data
    # so as to estimate variance of combined sample without 
    # the mean difference contributing to the variance estimate.
 
    # First we compute some summary statistics from the data
    ObsCovs <- cov(rbind(group1 - rep(c(1,2), each=n), group2))
    ObsMeansGroup1 <- c(mean(group1[,1]), mean(group1[,2]))
    ObsMeansGroup2 <- c(mean(group2[,1]), mean(group2[,2]))

    # Second we extract the parameter estimates and matrix algebra results from the model
    Sigma<-defmeansresult@matrices$S@values[1:2,1:2]
    Mu<-defmeansresult@matrices$M@values[1:2]
    beta<-defmeansresult@matrices$A@values[1:2,3]

    # Third, we check to see if things are more or less equal
    omxCheckCloseEnough(ObsCovs,Sigma,.01)
    omxCheckCloseEnough(ObsMeansGroup1,as.vector(Mu+beta),.001)
    omxCheckCloseEnough(ObsMeansGroup2,as.vector(Mu),.001)


