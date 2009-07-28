Regression, Path Specification
=====================================

Our next example will show how regression can be carried out from a path-centric structural modeling perspective. This example is in three parts; a simple regression, a multiple regression, and multivariate regression. There are two versions of each example are available; one with raw data, and one where the data is supplied as a covariance matrix and vector of means. These examples are availabe in the following files:

*SimpleRegression_PathCov.R
*SimpleRegression_PathRaw.R
*MultipleRegression_PathCov.R
*MultipleRegression_PathRaw.R
*MultivariateRegression_PathCov.R
*MultivariateRegression_PathRaw.R

Simple Regression
-----------------

We begin with a single dependent variable (y) and a single independent variable (x). The relationship between these variables takes the following form:

.. math::
   :nowrap:
\begin{eqnarray*} 
y &=& \beta_0 + \beta_1 * x + \epsilon
\end{eqnarray*}

In this model, the mean of y is dependent on both regression coefficients (and by extension, the mean of x). The variance of y depends on both the residual variance and the product of the regression slope and the variance of x. This model contains five parameters from a structural modeling perspective ($\beta_0$, $\beta_1$, $\sigma^2_{\epsilon}$, and the mean and variance of x). We are modeling a covariance matrix with three degrees of freedom (two variances and one variance) and a means vector with two degrees of freedom (two means). Because the model has as many parameters (5) as the data have degrees of freedom, this model is fully saturated.

Data
----

Our first step to running this model is to put include the data to be analyzed. The data must first be placed in a variable or object. For raw data, this can be done with the read.table function. The data provided has a header row, indicating the names of the variables.

.. code-block:: r
myRegDataRaw<-read.table("myRegData.txt",header=TRUE)

The names fo the variables provided by the header row can be displayed with the names() function.

.. code-block:: r
> names(myRegDataRaw)
[1] "w" "x" "y" "z"

As you can see, our data has four variables in it. However, our model only contains two variables, x and y. To use only them, we'll select only the variables we want and place them back into our data object. That can be done with the R code below.

.. We can refer to individual rows and columns of a data frame or matrix using square brackets, with selected rows referenced first and selected columns referenced second, separated by a comma. In the code below, we select all rows (there is no selection operator before the comma) and only columns x and y. As we are selecting multiple columns, we use the c() function to concatonate or connect those two names into one object.

.. code-block:: r
SimpleDataRaw <- myRegDataRaw[,c("x","y")]

For covariance data, we do something very similar. We create an object to house our data. Instead of reading in raw data from an external file, we can also include a covariance matrix. This requires the matrix() function, which needs to know what values are in the covariance matrix, how big it is, and what the row and column names are. As our model also references means, we'll include a vector of means in a separate object. Data is selected in the same way as before.

.. We'll select variables in much the same way as before, but we must now select both the rows and columns of the covariance matrix.  This means vector doesn't include names, so we'll just select the second and third elements of that vector.

.. code-block:: r
require(OpenMx)
myRegDataCov <- matrix(
    c(0.808,-0.110, 0.089, 0.361,
     -0.110, 1.116, 0.539, 0.289,
      0.089, 0.539, 0.933, 0.312,
      0.361, 0.289, 0.312, 0.836),
    nrow=4,
    dimnames=list(
        c("w","x","y","z"),
        c("w","x","y","z"))
)
 
SimpleDataCov <- myRegDataCov[c("x","y"),c("x","y")]	
 
myRegDataMeans <- c(2.582, 0.054, 2.574, 4.061)
 
SimpleDataMeans <- myRegDataMeans[c(2,3)]
	
Specifying the Model
--------------------

The following code contains all of the components of our model. All objects required for estimation (data, paths, and a model type) are included in their own arguments or functions. This code uses the ``mxModel`` function to create an ``MxModel`` object, which we'll then run.

.. code-block:: r
 uniRegModel <- mxModel("Simple Regression -- Path Specification", 
     type="RAM",
     mxData(
         data=myRegDataCov, 
         type="cov", 
         numObs=100,
         means=myRegDataMeans 
     ),
     manifestVars=c("x", "y"),
     # variances paths
     mxPath(
         from=c("x", "y"), 
         arrows=2,
         free=TRUE, 
         values = c(1, 1),
         labels=c("varx", "residual")
     ),
     # regression weights
     mxPath(
         from="x",
         to="y",
         arrows=1,
         free=TRUE,
         values=1,
         labels="beta1"
     ), 
     # means and intercepts
     mxPath(
         from="one",
         to=c("x", "y"),
         arrows=1,
         free=TRUE,
         values=c(1, 1),
         labels=c("meanx", "beta0")
     )
 ) # close model

