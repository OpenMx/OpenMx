Definition Variables, Matrix Specification
==========================================

This example will demonstrate the use of OpenMx definition variables with the implementation of a simple two group dataset.  What are definition variables?  Essentially, definition variables can be thought of as observed variables which are used to change the statistical model on an individual case basis.  In essence, it is as though one or more variables in the raw data vectors are used to specify the statistical model for that individual.  Many different types of statistical model can be specified in this fashion; some  are readily specified in standard fashion, and some that cannot.  To illustrate, we implement a two-group model.  The groups differ in their means but not in their variances and covariances.  This situation could easily be modeled in a regular multiple group fashion - it is only implemented using definition variables to illustrate their use.  The results are verified using summary statistics and an Mx 1.0 script for comparison is also available.

Mean Differences
----------------

The scripts are presented here

* DefinitionMeans_MatrixRaw.R
* DefinitionMeans_MatrixRaw.mx

Statistical Model
^^^^^^^^^^^^^^^^^

Algebraically, we are going to fit the following model to the observed x and y variables:

.. math::
   :nowrap:
   
   \begin{eqnarray*} 
   x_{i} = \mu_{x} + \beta_x * def + \epsilon_{xi}
   y_{i} = \mu_{y} + \beta_y * def + \epsilon_{yi}
   \end{eqnarray*}

where the residual sources of variance, :math:`\epsilon_{xi}` and :math:`\epsilon_{yi}` covary to the extent :math:`\rho`.  So, the task is to estimate: the two means :math:`\mu_{x}` and :math:`\mu_{y}`; the deviations from these means due to belonging to the group identified by having def set to 1 (as opposed to zero), :math:`\beta_{x}` and :math:`\beta_{y}`; and the parameters of the variance covariance matrix: cov(:math:`\epsilon_{x},\epsilon_{y}`) which we will call :math:`\Sigma` or simply "S" in the R script.

Data Simulation
^^^^^^^^^^^^^^^

Our first step to running this model is to simulate the data to be analyzed. Each individual is measured on two observed variables, x and y, and a third variable "def" which denotes their group membership with a 1 or a 0.  These values for group membership are not accidental, and must be adhered to in order to obtain readily interpretable results.  Other values such as 1 and 2 would yield the same model fit, but would make the interpretation more difficult.  

.. code-block:: r

    library(MASS) 

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

The objects *y* and *def* might be combined in a data frame.  However, in this case we won't bother to do it externally, and simply paste them together in the mxData function call.

Model Specification
^^^^^^^^^^^^^^^^^^^

The following code contains all of the components of our model. Before running a model, the OpenMx library must be loaded into R using either the ``require()`` or ``library()`` function. This code uses the ``mxModel`` function to create an ``mxModel`` object, which we'll then run.  Note that all the objects required for estimation (data, matrices, and an objective function) are declared within the ``mxModel`` function.  This type of code structure is recommended for OpenMx scripts generally.

.. code-block:: r

    defMeansModel <- mxModel("DefinitionMeans", 
        mxFIMLObjective(
            covariance="Sigma",
            means="Mu"
        ), 

The first argument in an ``mxModel`` function has a special function. If an object or variable containing an ``MxModel`` object is placed here, then ``mxModel`` adds to or removes pieces from that model. If a character string (as indicated by double quotes) is placed first, then that becomes the name of the model. Models may also be named by including a ``name`` argument. This model is named ``"DefinitionMeans"``.

The second argument in this mxModel call is itself a function. It declares that the objective function to be optimized is full information maximum likelihood (FIML) under normal theory, which is tagged as ``mxFIMLObjective``.  There are in turn two arguments to this function: the covariance matrix ``Sigma`` and the mean vector ``Mu``.  These matrices will be defined later in the mxModel function call.

Next, we declare where the data are, and their type, by creating an ``MxData`` object with the ``mxData`` function.
This piece of code creates an ``MxData`` object. It first references the object where our data are, then uses the ``type`` argument to specify that this is raw data. Analyses using definition variables have to use raw data, so that the model can be specified on an individual data vector level.

.. code-block:: r

        mxData((
            observed=data.frame(y,def)), 
            type="raw"
        ),

Model specification is carried out using ``mxMatrix`` functions to create matrices for the model. In the present case, we need four matrices.  First is the predicted covariance matrix, *Sigma*.  Next, we use three matrices to specify the model for the means.  First is *M* which corresponds to estimates of the means for individuals with definition variables with values of zero.  Individuals with definition variable values of 1 will have the value in *M* along with the value in the matrix *beta*.  So both matrices are of size 1x2 and both contain two free parameters.  There is a separate deviation for each of the variables, which will be estimated in the elements 1,1 and 1,2 of the *beta* matrix.  Last, but by no means least, is the matrix *def* which contains the definition variable.  The variable *def* in mxData data frame is referred to as ``data.def``.  In the present case, the definition variable contains a 1 for group 1, and a zero otherwise.  

.. code-block:: r

        mxMatrix(
            type="Symm", 
            nrow=2, 
            ncol=2, 
            free=TRUE, 
            values=c(1, 0, 1), 
            dimnames=list(selvars,selvars), 
            name="Sigma"
        ),
        mxMatrix(
            type="Full", 
            nrow = 2, 
            ncol = 1, 
            free=TRUE, 
            dimnames=list(selvars,NULL), 
            name = "M"
        ),
        mxMatrix(
            type="Full", 
            nrow=2, 
            ncol=1, 
            free=TRUE, 
            values=c(0, 0),
            dimnames=list(selvars,NULL), 
            name="beta"
        ),
        mxMatrix(
            type="Full", 
            nrow=2, 
            ncol=1, 
            free=FALSE, 
            labels=c("data.def"),
            dimnames=list(selvars,NULL), 
            name="def"
        ),

The trick - commonly used in regression models - is to multiply the *beta* matrix by the *def* matrix.  This multiplication is effected using an mxAlgebra function call:

.. code-block:: r

        mxAlgebra(
            expression=M+beta*def, 
            name="Mu", 
            dimnames=list(selvars,NULL)
        )
    )

The result of this algebra is named *Mu``, and this handle is referred to in the mxFIMLObjective function call.  We can then run the model and examine the output with a few simple commands.

Model Fitting
^^^^^^^^^^^^^^

.. code-block:: r

    # Run the model
    defMeansFit <- mxRun(defMeansModel)
    defMeansFit@matrices
    defMeansFit@algebras

It is possible to compare the estimates from this model to some summary statistics computed from the data:

.. code-block:: r

    # Compare OpenMx estimates to summary statistics computed from raw data.
    # Note that to calculate the common variance, 
    # group 1 has the 1 and 2 subtracted from every Xi and Yi in the sample
    # data, so as to estimate variance of combined sample without the mean correction.
 
    # First we compute some summary statistics from the data
    ObsCovs<-cov(rbind(group1-rep(c(1,2),each=n),group2))
    ObsMeansGroup1<-c(mean(group1[,1],mean(group1[,2]))
    ObsMeansGroup2<-c(mean(group2[,1],mean(group2[,2]))
 
    # Second we extract the parameter estimates and matrix algebra results from the model
    Sigma<-run@matrices$Sigma@values
    Mu<-run@algebras$Mu@result
    M<-run@matrices$M@values
    beta<-run@matrices$beta@values
 
    # Third, we check to see if things are more or less equal
    omxCheckCloseEnough(ObsCovs,Sigma,.01)
    omxCheckCloseEnough(ObsMeansGroup1,as.vector(M+beta),.001)
    omxCheckCloseEnough(ObsMeansGroup2,as.vector(Mu),.001)


