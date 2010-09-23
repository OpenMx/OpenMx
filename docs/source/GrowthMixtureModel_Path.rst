
Growth Mixture Modeling, Path Specification
===========================================

This example will demonstrate how to specify a growth mixture model using path specification. Unlike other examples, this application will not be demonstrated with covariance data, as this model can only be fit to raw data. The script for this example can be found in the following file:

** http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/GrowthMixtureModel_Path.R

A parallel example using matrix specification can be found here:

** http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/GrowthMixtureModel_Matrix.R

The latent growth curve used in this example is the same one fit in the latent growth curve example. Path and matrix versions of that example for raw data can be found here: 

** http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/LatentGrowthCurveModel_PathRaw.R

** http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/LatentGrowthCurveModel_MatrixRaw.R

Mixture Modeling
----------------

Mixture modeling is an approach where data are assumed to be governed by some type of mixture distribution. This includes a large class of models, including many varieties of mixture modeling, latent class analysis and related models with binary or categorical latent variables. This example will demonstrate a growth mixture model, where change over time is modeled with a linear growth curve and the distribution of latent intercepts and slopes is governed by a mixture of two distributions. The model can thus be described as a combination of two growth curves, weighted by a class proportion variable, as shown below.

.. math::
   :nowrap:
   
   \begin{eqnarray*} 
   x_{ij} = p_1 (Intercept_{i1} + \lambda_1 Slope_{i1} + \epsilon) + p_2 (Intercept_{i2} + \lambda_2 Slope_{i2} + \epsilon)
   \end{eqnarray*}

To scale the class proportion variable as a probability, it must be scaled such that it is strictly positive and the set of all class probabilities sum to a value of one.

.. math::
   :nowrap:

   \begin{eqnarray*} 
   \sum_{i=1}^k p_i = 1 
   \end{eqnarray*}

Data
^^^^
The data for this example can be found in the data object ``myGrowthMixtureData``. These data contain five time ordered variables named ``x1`` through ``x5``, just like the growth curve demo mentioned previously. It is important to note that raw data is required for mixture modeling, as moment matrices do not contain all of the information required to estimate the model. 

.. code-block:: r

	data(myGrowthMixtureData)
	names(myGrowthMixtureData)

Model Specification
^^^^^^^^^^^^^^^^^^^

Specifying a mixture model can be categorized into two general phases. The first phase of model specification pertains to creating the models for each class. The second phase specifies the way those classes are mixed. In OpenMx, this is done using a model tree. Each class is created as a separate ``MxModel`` object, and those class-specific models are all placed into a larger or parent model. The parent model contains the class proportion parameter(s) and the data. 

Creating the class-specific models is done the same way as every other model. We'll begin by specifying the model for the first class using the ``mxPath`` function. The code below specifies a five-occasion linear growth curve, virtually identical to the one in the linear growth curve example referenced above. The only changes made to this model are the names of the free parameters; the means, variances and covariance of the intercept and slope terms are now followed by the number 1 to distinguish them from free parameters in the other class.

.. code-block:: r

	class1 <- mxModel("Class1", 
	    type="RAM",
	    manifestVars=c("x1","x2","x3","x4","x5"),
	    latentVars=c("intercept","slope"),
		# residual variances
	    mxPath(
	    	from=c("x1","x2","x3","x4","x5"), 
	        arrows=2,
	        free=TRUE, 
	        values = c(1, 1, 1, 1, 1),
	        labels=c("residual","residual","residual","residual","residual")
	    ),
  	  # latent variances and covariance
	    mxPath(
	    	from=c("intercept","slope"), 
	        arrows=2,
	        all=TRUE,
	        free=TRUE, 
	        values=c(1, .4, .4, 1),
	        labels=c("vari1", "cov1", "cov1", "vars1")
	    ),
	    # intercept loadings
	    mxPath(
	    	from="intercept",
	        to=c("x1","x2","x3","x4","x5"),
	        arrows=1,
	        free=FALSE,
	        values=c(1, 1, 1, 1, 1)
	    ),
	    # slope loadings
	    mxPath(
	    	from="slope",
	        to=c("x1","x2","x3","x4","x5"),
	        arrows=1,
	        free=FALSE,
	        values=c(0, 1, 2, 3, 4)
	    ),
	    # manifest means
	    mxPath(from="one",
	        to=c("x1", "x2", "x3", "x4", "x5"),
	        arrows=1,
	        free=FALSE,
	        values=c(0, 0, 0, 0, 0)
	    ),
	    # latent means
	    mxPath(from="one",
	        to=c("intercept", "slope"),
	        arrows=1,
	        free=TRUE,
	        values=c(0, -1),
	        labels=c("meani1", "means1")
	    )
	) # close model
	
We could create the model for our second class by copy and pasting the code above, but that can yield needlessly long scripts. We can also use the ``mxModel`` function to edit an existing model object, allowing us to change only the parameters that vary across classes. The ``mxModel`` call below begins with an existing ``MxModel`` object (``class1``) rather than a model name. The subsequent ``mxPath`` functions add new paths to the model, replacing any existing paths that describe the same relationship. As we did not give the model a name at the beginning of the ``mxModel`` function, we must use the ``name`` argument to identify this model by name.

.. code-block:: r

	class2 <- mxModel(class1,
		# latent variances and covariance
	    mxPath(
	    	from=c("intercept","slope"), 
	        arrows=2,
	        all=TRUE,
	        free=TRUE, 
	        values=c(1, .5, .5, 1),
	        labels=c("vari2", "cov2", "cov2", "vars2")
	    ),
	    # latent means
	    mxPath(from="one",
	        to=c("intercept", "slope"),
	        arrows=1,
	        free=TRUE,
	        values=c(5, 1),
	        labels=c("meani2", "means2")
	    ),
		name="Class2"
	) # close model

We must make one other change to our class-specific models before creating the parent model that will contain them. The objective function for each of the class-specific models must return the likelihoods for each individual rather than the default log likelihood for the entire sample. OpenMx objective functions that handle raw data have the option to return a vector of likelihoods for each row rather than a single likelihood value for the dataset. This option can be accessed either as an argument in a function like ``mxRAMObjective`` or ``mxFIMLObjective`` or with the syntax below.

.. code-block:: r

	class1@objective@vector <- TRUE
	class2@objective@vector <- TRUE
	
While the class-specific models can be specified using either path or matrix specification, the class proportion parameter must be specified using a matrix, though it can be specified a number of different ways. The code below demonstrates one method of specifying class proportion parameters as probabilities. 

The matrix in the object ``classP`` contains two elements representing the proportion of the sample in each of the two classes, while the object ``classA`` contains an ``MxAlgebra`` that scales this proportion as a probability. Placing bounds on the class probabilities matrix constrains each of the probabilities to be between zero and one, while the algebra defines the probability of being in class 2 to be 1 minus the probability of being in class 1. This ensures that the sum of the class probabilities is 1. Notice that the second element of the class probability matrix is constrained to be equal to the result of the ``mxAlgebra`` statement. The brackets in the ``mxMatrix`` function are required; the second element in the "classProbs" object is actually constrained to be equal to the first row and first column of the ``MxAlgebra`` object "pclass2", which evaluates to a 1 x 1 matrix.

.. code-block:: r

	classP <- mxMatrix("Full", 2, 1, free=c(TRUE, FALSE), 
	          values=.2, lbound=0.001, ubound=0.999,
	          labels = c("pclass1", "pclass2[1,1]"), name="classProbs")

	classA <- mxAlgebra(1-pclass1, name="pclass2")
	
The above code creates one free parameter for class probability ("pclass1") and one fixed parameter, which is the result of an algebra ("pclass2"). There are at least two other ways to specify this class proportion parameter, each with benefits and drawbacks. One could create two free parameters named "pclass1" and "pclass2" and constrain them using the ``mxConstraint`` function. This approach is relatively straightforward, but comes at the expense of standard errors. Alternatively, one could omit the algebra and fix "pclass2" to a specific value. This would make model specification easier, but the resulting "pclass1" parameter would not be scaled as a probability.

Finally, we can specify the mixture model. We must first specify the model's -2 log likelihood function defined as:

.. math::
   :nowrap:
   
   \begin{eqnarray*} 
   -2LL = -2 * \sum \log (p_1 l_{1i} + p_2 l_{2i})
   \end{eqnarray*}
	
This is specified using an ``mxAlgebra`` function, and used as the argument to the ``mxAlgebraObjective`` function. Then the objective function, matrices and algebras used to define the mixture distribution, the models for the respective classes and the data are all placed in one final ``mxModel`` object, shown below.	

.. code-block:: r

	algObj <- mxAlgebra(-2*sum(
	          log(pclass1%x%Class1.objective + pclass2%x%Class2.objective)), 
	          name="mixtureObj")

	obj <- mxAlgebraObjective("mixtureObj")
	
	gmm <- mxModel("Growth Mixture Model",
		mxData(
	    	observed=myGrowthMixtureData,
	        type="raw"
	    ),
	    class1, class2,
	    classP, classA,
	    algObj, obj
		)      

	gmmFit <- mxRun(gmm)

	summary(gmmFit)

Multiple Runs
^^^^^^^^^^^^^

The results of a mixture model can sometimes depend on starting values. It is a good idea to run a mixture model with a variety of starting values to make sure results you find are not the result of a local minimum in the likelihood space.

One way to access the starting values in a model is by using the ``omxGetParameters`` function. This function takes an existing model as an argument and returns the names and values of all free parameters. Using this function on our growth mixture model, which is stored in an objected called ``gmm``, gives us back the starting values we specified above.

.. code-block:: r

        omxGetParameters(gmm)
            pclass1 residual    vari1     cov1    vars1   meani1   means1    vari2     cov2    vars2   meani2 
            	0.2      1.0      1.0      0.4      1.0      0.0     -1.0      1.0      0.5      1.0      5.0 
            means2 
            	1.0

A companion function to ``omxGetParameters`` is ``omxSetParameters``, which can be used to alter one or more named parameters in a model. This function can be used to change the values, freedom and labels of any parameters in a model, returning an MxModel object with the specified changes. The code below shows how to change the residual variance starting value from 1.0 to 0.5. Note that the output of the ``omxSetParameters`` function is placed back into the object ``gmm``.

.. code-block:: r

		gmm <- omxSetParameters(gmm, labels="residual", values=0.5)

The MxModel in the object ``gmm`` can now be run and the results compared with other sets of staring values. Starting values can also be sampled from distributions, allowing users to automate starting value generation, which is demonstrated below. The ``omxGetParameters`` function is used to find the names of the free parameters and define three matrices: a matrix ``input`` that holds the starting values for any run; a matrix ``output`` that holds the converged values of each parameter; and a matrix ``fit`` that contains the -2 log likelihoods and other relevant model fit statistics. Each of these matrices contains one row for every set of starting values. A ``for`` loop repeatedly generates starting values (from a set of uniform distributions using ``runif``), runs the model with those starting values and places the starting values, final estimates and fit statistics in the ``input``, ``output`` and ``fit`` matrices, respectively.

.. code-block:: r

	trials <- 20

	omxGetParameters(gmm)

	parNames <- names(omxGetParameters(gmm))
	
	input <- matrix(NA, trials, length(parNames))
	dimnames(input) <- list(c(1: trials), c(parNames))

	output <- matrix(NA, trials, length(parNames))
	dimnames(output) <- list(c(1: trials), c(parNames))

	fit <- matrix(NA, trials, 4)
	dimnames(fit) <- list(c(1: trials), c("Minus2LL", "Status", "Iterations", "pclass1"))
	
	for (i in 1: trials){
		cp <- runif(1, 0.1, 0.9) # class probability
		v  <- runif(5, 0.1, 5.0) # variance terms
		cv <- runif(2,-0.9, 0.9) # covariances (as correlations)
		m  <- runif(4,-5.0, 5.0) # means
		cv <- cv*c(sqrt(v[2]*v[3]), sqrt(v[4]*v[5])) #rescale covariances
	
		temp1 <- omxSetParameters(gmm,
			labels=parNames,
			values=c(
				cp, # class probability
				v[1], 
				v[2], cv[1], v[3], m[1], m[2],
				v[4], cv[2], v[5], m[3], m[4]
				)
				)
		
				temp1@name <- paste("Starting Values Set", i)
		
				temp2 <- mxRun(temp1, unsafe=TRUE)
	
				input[i,] <- omxGetParameters(temp1)
				output[i,] <- omxGetParameters(temp2)
				fit[i,] <- c(
				temp2@output$Minus2LogLikelihood,
				temp2@output$status[[1]],
				temp2@output$iterations,
				temp2@output$estimate[1]
				)
			}
	
Viewing the contents of the ``fit`` matrix shows the -2 log likelihoods for each of the runs, as well as the convergence status, number of iterations and class probabilities, shown below.

.. code-block:: r

	fit
	   Minus2LL Status Iterations   pclass1
	1  8739.050      0         41 0.3991078
	2  8739.050      0         40 0.6008913
	3  8739.050      0         44 0.3991078
	4  8739.050      1         31 0.3991079
	5  8739.050      0         32 0.3991082
	6  8739.050      1         34 0.3991089
	7  8966.628      0         22 0.9990000
	8  8966.628      0         24 0.9990000
	9  8966.628      0         23 0.0010000
	10 8966.628      1         36 0.0010000
	11 8963.437      6         25 0.9990000
	12 8966.628      0         28 0.9990000
	13 8739.050      1         47 0.6008916
	14 8739.050      1         36 0.3991082
	15 8739.050      0         43 0.3991076
	16 8739.050      0         46 0.6008948
	17 8739.050      1         50 0.3991092
	18 8945.756      6         50 0.9902127
	19 8739.050      0         53 0.3991085
	20 8966.628      0         23 0.9990000

There are several things to note about the above results. First, the minimum -2 log likelihood was reached in 12 of 20 sets of staring values, all with NPSOL statuses of either zero (seven times) or one (five times). Additionally, the class probabilities are equivalent within five digits of precision, keeping in mind that no the model as specified contains no restriction as to which class is labeled "class 1" (probability equals .3991) and "class 2" (probability equals .6009). The other eight sets of starting values showed higher -2 log likelihood values and class probabilities at the set upper or lower bounds, indicating a local minimum. We can also view this information using R's ``table`` function.

.. code-block:: r

	table(round(fit[,1], 3), fit[,2])
          
	           0 1 6
	  8739.05  7 5 0
	  8945.756 0 0 1
	  8963.437 0 0 1
	  8966.628 5 1 0

We should have a great deal of confidence that the solution with class probabilities of .399 and .601 is the correct one.

Multicore Estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OpenMx supports multicore processing through the ``snowfall`` library, which is described in the "Multicore Execution" section of the documentation and in the following demo:

** http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/BootstrapParallel.R

Using multiple processors can greatly improve processing time for model estimation when a model contains independent submodels. While the growth mixture model in this example does contain submodels (i.e., the class specific models), they are not independent, as they both depend on a set of shared parameters ("residual", "pclass1").

However, multicore estimation can be used instead of the ``for`` loop in the above section for testing alternative sets of starting values. Instead of changing the starting values in the ``gmm`` object repeatedly, multiple copies of the model contained in ``gmm`` must be placed into parent or container model. Either the above ``for`` loop or a set of "apply" statements can be used to generate the model.
