Beginners Guide to OpenMx
=========================

This document is a basic introduction to OpenMx.  It assumes that the reader has installed the R statistical programming language [http://www.r-project.org/] and the OpenMx library for R [http://openmx.psyc.virginia.edu].  Detailed introductions to R can be found on the internet.  We recommend [http://faculty.washington.edu/tlumley/Rcourse] for a short course but Google search for 'introduction to R' provides many options.

The OpenMx scripts for the examples in this guide are available in the following files:

* http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/OneFactorPathDemo.R
* http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/OneFactorMatrixDemo.R
* http://www.vipbg.vcu.edu/OpenMx/html/NewBeginnersGuide.R

.. _BasicIntroduction:

Basic Introduction 
------------------

The main function of OpenMx is to design a statistical model which can be used to test a hypothesis with a set of data.  The way to do this is to use one of the functions written in the R language as part of the OpenMx package.  The function to create an Mx model is (you guessed it) ``mxModel()``.  Note that i) R is case sensitive and ii) the OpenMx package must be loaded before any of the OpenMx functions are used (this only has to be done once in an R session).

..
   DO NOT EXECUTE

.. cssclass:: input
   
   .. code-block:: r
       
        require(OpenMx)

We start by building a model with the ``mxModel()`` function, which takes a number of arguments.  We explain each of these below.  Usually we save the result of applying the ``mxModel()`` function as an R object, here called *myModel*.  

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r
       
        myModel <- mxModel() 

This R object has a special class named MxModel. We might read the above line 'myModel gets mxModel left paren right paren'. In this case we have provided no arguments in the parentheses, but R has still created an empty MxModel object *myModel*. Obviously this model needs arguments to do anything useful, but just to get the idea of the process of 'Model Building - Model Fitting' here is how we would go about fitting this model.  The *myModel* object becomes the argument of the OpenMx function ``mxRun()`` to fit the model.  The result is stored in a new R object, *myModelRun* of the same class as the previous object but updated with the results from the model fitting.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r
       
        myModelRun <- mxRun(myModel) 

A model typically contains data, free and/or fixed parameters and some function to be applied to the data according to the model.  Models can be build in a variety of ways, due to the flexibility of OpenMx which it inherited from classic **Mx** - for those of you who familiar with that software - and further expanded.

One possibility for building models is to create all the elements as separate R objects, which are then combined by calling them up in the ``mxModel()`` statement.  We will refer to this approach as the **piecewise style**.

A second type is to start by creating an ``mxModel()`` with one element, and then taking this model as the basis for the next model where you add another element to it and so on.  This also provides a good way to make sure that each element is syntactically correct before adding anything new.  We will call this approach the **iterative/recursive/stepwise style**.

Another approach more closely resembles the traditional classic **Mx** approach, where you specify all the elements of the model at once, all as arguments of the ``mxModel()``, separated by comma's.  While this approach is often more compact and works well for scripts that run successfully, it is not the most recommended approach for debugging a new script.  We will refer to this approach as the **classic style**.


A First mxModel
----------------

To introduce the model fitting process in OpenMx, we will present the basics of several OpenMx functions which can be used to write a simple model and view its contents.

Matrix Creation
^^^^^^^^^^^^^^^

Although ``mxModel()`` can have a range of arguments, we will start with the most simple one.  Models are fitted to data which must be in numeric format (for continuous data) or factor format (for ordinal data).  Here we consider continuous data.  Numbers (data/parameter estimates) are typically put into matrices, except for fixed constants.  The function created to put numbers into matrices is (unsurprisingly) ``mxMatrix()``.  Here we start with a basic matrix call and make use of only some of its possible arguments. All arguments are separated by comma's. To make it clear and explicit, we will include the names of the arguments, although that is optional if the arguments are included in the default order.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r
       
        myAmatrix <- mxMatrix(type="Full", nrow=1, ncol=1, values=4, name="Amatrix")
    
The above call to the ``mxMatrix()`` function has five arguments.  The ``type`` and ``name`` arguments are alphanumeric and therefore their values are in quotes.  The ``nrows``, ``ncols`` and ``values`` arguments are numeric, and refer respectively to the number of rows, the number of columns of the matrix and the value for the (in this case only one) element of the matrix.

Matrix Contents
^^^^^^^^^^^^^^^

Once you have run/executed this statement in R, a new R object has been created, namely *myAmatrix*.  When you view its contents, you'll notice it has a special class of object, made by OpenMx, called an MxMatrix object.  This object has a number of attributes, all of which are listed when you call up the object.  

    ..  code-block:: r

        > myAmatrix
        FullMatrix 'Amatrix' 
        
        $labels: No labels assigned.
        
        $values
          [,1]
        [1,]    4
        
        $free: No free parameters.
        
        $lbound: No lower bounds assigned.
        
        $ubound: No upper bounds assigned.

Most of these attributes start with the ``$`` symbol.  The contents of a particular attribute can be displayed by typing the name of the R object followed by the ``$`` symbol and the name of the attribute, for example here we're displaying the values of the matrix *myAmatrix*
   
   .. code-block:: r
   
        > myAmatrix$values
               [,1]
          [1,]    4

Note that the attribute ``name`` is part of the header of the output but is not displayed as an ``$`` attribute.  However, it does exist as one and can be seen by typing
   
   .. code-block:: r
   
        > myAmatrix$name
        [1] "Amatrix"

Wait a minute, this is confusing.  The matrix has a name, here "Amatrix", and the R object to represent the matrix has a name, here "myAmatrix".  Remember that when you call up *myAmatrix* you get the contents of the entire MxMatrix R object.  When you call up "Amatrix", you get 

    .. code-block:: r

        Error: object 'Amatrix' not found   

unless you had previously created another R object with that same name.  Why do we need two names?  The matrix name (here, "Amatrix") is used within OpenMx when performing an operation on this matrix using algebra (see below) or manipulating/using the matrix in any way within a model.  When you want to manipulate/use/view the matrix outside of OpenMx, or build a model by building each of the elements as R objects in the 'piecewise' approach, you use the R object name (here, *myAmatrix*).  Let's clarify this with an example.  

Model Creation
^^^^^^^^^^^^^^

First, we will build a model *myModel1* with just one matrix.  Obviously that is not very useful but it does serve to introduce the sequence of creating a model and running it.  

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        myModel1     <- mxModel( mxMatrix(type="Full", nrow=1, ncol=1, values=4, name="Amatrix") )
                    
Model Execution
^^^^^^^^^^^^^^^^

The ``mxRun()`` function will run a model through the optimizer.  The return value of this function is an identical MxModel object, with all the free parameters - in case there are any - in the elements of the matrices of the model assigned to their final values.                    
                    
.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r
   
        myModel1Run  <- mxRun(myModel1)

Model Contents
^^^^^^^^^^^^^^

Note that we have saved the result of applying ``mxRun()`` to *myModel1* into a new R object, called *myModel1Run* which is of the same class as *myModel1* but with values updated after fitting the model.  Note that the MxModel is automatically given a name 'untitled2' as we did not specify a ``name`` argument for the ``mxModel()`` function.
   
   .. code-block:: r
    
        >     myModel1Run
        MxModel 'untitled2' 
        type : default 
        $matrices : 'Amatrix' 
        $algebras :  
        $constraints :  
        $intervals :  
        $latentVars : none
        $manifestVars : none
        $data : NULL
        $submodels :  
        $expectation : NULL 
        $fitfunction : NULL 
        $compute : NULL 
        $independent : FALSE 
        $options :  
        $output : TRUE

As you can see from viewing the contents of the new object, the current model only uses two of the arguments, namely ``$matrices`` and ``$output``.  Given the matrix was specified within the mxModel, we can explore its arguments by extending the level of detail as follows.

   .. code-block:: r

        > myModel1Run$matrices
        $Amatrix
        FullMatrix 'Amatrix' 
    
        $labels: No labels assigned.
    
        $values
             [,1]
        [1,]    4
    
        $free: No free parameters.
    
        $lbound: No lower bounds assigned.
    
        $ubound: No upper bounds assigned.
    
This lists all the matrices within the MxModel *myModel1Run*.  In the current case there is only one.  If we want to display just a specific argument of that matrix, we first add a dollar sign ``$``, followed by the name of the matrix, and an ``$`` sign prior to the required argument.  Thus both arguments within an object and specific elements of the same argument type are preceded by the ``$`` symbol.

    .. code-block:: r

        > myModel1run$matrices$Amatrix$values
              [,1]
         [1,]    4

It is also possible to omit the ``$matrices`` part and use the more succinct ``myModel1Run$Amatrix$values``.

Similarly, we can inspect the output which also includes the matrices in ``$matrices``, but only displays the values.  Furthermore, the output will list algebras (``$algebras``), model expectations (``$expectations``), status of optimization (``$status``), number of evaluations (``$evaluations``), openmx version (``$mxVersion``), and a series of time measures of which the CPU time might be most useful (``$cpuTime``).

    .. code-block:: r

        > myModel1Run$output
        $matrices
        $matrices$untitled2.Amatrix
             [,1]
        [1,]    4

        ....
        $mxVersion
        [1] "999.0.0-3297"

        $frontendTime
        Time difference of 0.05656791 secs

        $backendTime
        Time difference of 0.003615141 secs

        $independentTime
        Time difference of 3.385544e-05 secs

        $wallTime
        Time difference of 0.0602169 secs

        $timestamp
        [1] "2014-04-10 09:53:37 EDT"

        $cpuTime
        Time difference of 0.0602169 secs

Alternative 
^^^^^^^^^^^

Now let's go back to the model *myModel1* for a minute.  We specified the matrix "Amatrix" within the model.  Given we had previously saved the "Amatrix" in the *myAmatrix* object, we could have just used the R object as the argument of the model as follows.  Here we're adding one additional element to the ``MxModel()`` object, namely the ``name`` argument

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        myModel2     <- mxModel(myAmatrix, name="model2")
        myModel2Run  <- mxRun(myModel2)

You can verify for yourself that the contents of *myModel2* is identical to that of *myModel1*, and the same applies to *myModel1Run* and *myModel2Run*, and as a result to the matrix contained in the model.  The value of the matrix element is still 4, both in the original model and the fitted model, as we did not manipulate the matrix in any way.  We refer to this alternative style of coding as **iterative**.

Algebra Creation
^^^^^^^^^^^^^^^^

Now, let's take it one step further and use OpenMx to evaluate some matrix algebra.  It will come as a bit of a shock to learn that the OpenMx function to specify an algebra is called ``mxAlgebra()``.  Its main argument is the ``expression``, in other words the matrix algebra formula you want to evaluate.  In this case, we're simply adding 1 to the value of the matrix element, providing a name for the matrix "Bmatrix" and then save the new matrix as *myBmatrix*.  Note that the matrix we are manipulating is the "Amatrix", the name given to the matrix within OpenMx.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        myBmatrix    <- mxAlgebra(expression=Amatrix+1, name="Bmatrix")
    
Algebra Contents
^^^^^^^^^^^^^^^^

We can view the contents of this new matrix. Notice that the result has not yet computed, as we have not run the model yet.
    
   .. code-block:: r
    
        > myBmatrix
        mxAlgebra 'Bmatrix' 
        $formula:  Amatrix + 1 
        $result: (not yet computed) <0 x 0 matrix>
        dimnames: NULL

Built Model
^^^^^^^^^^^

Now we can combine the two statements - one defining the matrix, and the other defining the algebra - in one model, simply by separating them by a comma, and run it to see the result of the operation.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        myModel3     <- mxModel(myAmatrix, myBmatrix, name="model3")
        myModel3Run  <- mxRun(myModel3)

First of all, let us view *myModel3* and more specifically the values of the matrices within that model.  Note that the ``$matrices`` lists one matrix, "Amatrix", and that the ``$algebras`` lists another, "Bmatrix".  To view values of matrices created with the ``mxMatrix()`` function, the argument is ``$values``; for matrices created with the ``mxAlgebra()`` function, the argument is ``$result``.  Note that when viewing a specific matrix, you can omit the ``$matrices`` or the ``$algebras`` arguments.

   .. code-block:: r

        >     myModel3
        MxModel 'model3' 
        type : default 
        $matrices : 'Amatrix' 
        $algebras : 'Bmatrix' 
        $constraints :  
        $intervals :  
        $latentVars : none
        $manifestVars : none
        $data : NULL
        $submodels :  
        $expectation : NULL 
        $fitfunction : NULL 
        $compute : NULL 
        $independent : FALSE 
        $options :  
        $output : FALSE 

   .. code-block:: r

        > myModel3$Amatrix$values
             [,1]
        [1,]    4

   .. code-block:: r

        > myModel3$Bmatrix$result
        <0 x 0 matrix>

Fitted Model
^^^^^^^^^^^^

Given we're looking at the model *myModel3* before it is run, results of algebra have not been computed yet.  Let us see how things change after running the model and viewing *myModel3Run*.

   .. code-block:: r

        >     myModel3Run
        MxModel 'model3' 
        type : default 
        $matrices : 'Amatrix' 
        $algebras : 'Bmatrix' 
        $constraints :  
        $intervals :  
        $latentVars : none
        $manifestVars : none
        $data : NULL
        $submodels :  
        $expectation : NULL 
        $fitfunction : NULL 
        $compute : NULL 
        $independent : FALSE 
        $options :  
        $output : TRUE

   .. code-block:: r

        > myModel3Run$Amatrix$values
             [,1]
        [1,]    4
   
   .. code-block:: r

        > myModel3Run$Bmatrix$result
             [,1]
        [1,]    5

You will notice that the structure of the MxModel objects is identical, the value of the "Amatrix" has not changed, as it was a fixed element.  However, the value of the "Bmatrix" is now the result of the operation on the "Amatrix".  Note that we're here looking at the "Bmatrix" within the MxModel object *myModel3Run*.   Please verify that the original MxAlgebra objects *myBmatrix* and *myAmatrix* remain unchanged.  The ``mxModel()`` function call has made its own internal copies of these objects, and it is only these internal copies that are being manipulated.  In computer science terms, this is referred to as *pass by value*.


Pass By Value
^^^^^^^^^^^^^

Let us insert a mini-lecture on the R programming language.  Our experience has found that this exercise will greatly increase your understanding of the OpenMx language. 

As this is such a crucial concept in R (unlike many other programming languages), let us look at it in a simple R example.  We will start by assigning the value 4 to the object *avariable*, and then display it.  If we then add 1 to this object, and display it again, notice that the value of *avariable* has not changed.

.. cssclass:: input
..

   R Code
   
   .. code-block:: r

        > avariable <- 4
        > avariable
        [1] 4
        > avariable +1
        [1] 5
        > avariable
        [1] 4
    
Now we introduce a function, as OpenMx is a collection of purposely built functions.  The function takes a single argument (the object *number*), adds one to the argument *number* and assigns the result to *number*, and then returns the incremented number back to the user.  This function is given the name ``addone()``.  We then apply the function to the object *avariable*, as well as display *avariable*.  Thus, the objects *addone* and *avariable* are defined. The object assigned to *addone* is a function, while the value assigned to *avariable* is the number 4. 

.. cssclass:: input
..

   R Code
   
   .. code-block:: r

        > addone <- function(number) {
            number <- number + 1
            return(number)
            }

        > addone(avariable)
        [1] 5
        > avariable
        [1] 4

Note that it may be prudent to use the ``print()`` function to display the results back to the user.  When R is run from a script rather than interactively, results will not be displayed unless the function ``print()`` is used as shown below.

.. cssclass:: input
..

   R Code
   
   .. code-block:: r

        > print(addone(avariable))
        [1] 5
        > print(avariable)
        [1] 4

What is the result of executing this code? Try it. The correct results are 5 and 4.  But why is the object *avariable* still 4, even after the ``addone()`` function was called? The answer to this question is that R uses pass-by-value function call semantics.

In order to understand pass-by-value semantics, we must understand the difference between *objects* and *values*. The *objects* declared in this example are *addone*, *avariable*, and *number*.  The *values* refer to the things that are stored by the *objects*.  In programming languages that use pass-by-value semantics, at the beginning of a function call it is the *values* of the argument list that are passed to the function.  

The object *avariable* cannot be modified by the function ``addone()``.  If I wanted to update the value stored in the object, I would have needed to replace the expression as follows:

.. cssclass:: input
..

   R Code
   
   .. code-block:: r

        > print(avariable <- addone(avariable))
        [1] 5
        > print(avariable)
        [1] 5
    
Try it.  The updated example prints out 5 and 5.  The lesson from this exercise is that the only way to update a object in a function call is to capture the result of the function call [#f1]_.  This lesson is sooo important that we'll repeat it:

*the only way to update an object in a function call is to capture the result of the function call.*

R has several built-in types of values that you are familiar with: numerics, integers, booleans, characters, lists, vectors, and matrices. In addition, R supports S4 object values to facilitate object-oriented programming.  Most of the functions in the OpenMx library return S4 object values.  You must always remember that R does not discriminate between built-in types and S4 object types in its call semantics.  Both built-in types and S4 object types are passed by value in R (unlike many other languages).

.. rubric:: Footnotes

.. [#f1] There are a few exceptions to this rule, but you can be assured such trickery is not used in the OpenMx library.


Styles
------

In the beginning of the introduction, we discussed three styles of writing OpenMx code: the piecewise, stepwise and classic styles.  Let's take the most recent model and show how it can be written in these three styles.  

Piecewise Style
^^^^^^^^^^^^^^^

The style we used in *myModel3* is the piecewise style.  We repeat the different statements here for clarity

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        myAmatrix    <- mxMatrix(type="Full", nrow=1, ncol=1, values=4, name="Amatrix")
        myBmatrix    <- mxAlgebra(expression=Amatrix+1, name="Bmatrix")

        myModel3     <- mxModel(myAmatrix, myBmatrix, name="model3")
        myModel3Run  <- mxRun(myModel3)

Each argument of the ``mxModel()`` statement is defined separately first as independent R objects which are then combined in one model statement.

Stepwise Style
^^^^^^^^^^^^^^^

For the stepwise style, we start with an ``mxModel()`` with just one argument, as we originally did with the "Amatrix" in *myModel1*, as repeated below.  We could run this model to make sure it's syntactically correct.

.. cssclass:: input
..

    OpenMx Code

    .. code-block:: r

        myModel1     <- mxModel( mxMatrix(type="Full", nrow=1, ncol=1, values=4, name="Amatrix") )
        myModel1Run  <- mxRun(myModel1)
 
Then we would build a new model starting from the first model.  To do this, we invoke a special feature of the first argument of an ``mxModel()``.  If it is the name of a saved MxModel object, for example *myModel1*, the arguments of that model would be automatically included in the new model.  These arguments can be changed (or not) and new arguments can be added.  Thus, in our example, where we want to keep the "Amatrix" and add the "Bmatrix", our second model would look like this.  

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        myModel4     <- mxModel(myModel1,
                        mxAlgebra(expression=Amatrix+1, name="Bmatrix"),
                        name="model4"
                        )
        myModel4Run  <- mxRun(myModel4)
    
Note that we call it "model4", by adding a ``name`` argument to the ``mxModel()`` as to not overwrite our previous "model1".

Classic Style
^^^^^^^^^^^^^

The final style may be reminiscent of classic Mx.  Here we build all the arguments explicitly within one ``mxModel()``.  As a result only one R object is created prior to ``mxRun()`` ning the model.  This style is more compact than the others but harder to debug.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        myModel5     <- mxModel(
                        mxMatrix(type="Full", nrow=1, ncol=1, values=4, name="Amatrix"),
                        mxAlgebra(expression=Amatrix+1, name="Bmatrix"), 
                        name="model5"
                        )
        myModel5Run  <- mxRun(myModel5)

You may have seen an alternative version with the first argument in quotes.  In that case, that argument refers to the name of the model and not to a previously defined model.  Thus, the following specification is identical to the previous one.  Note also that it is not necessary to add the 'names' of the arguments, as long as the arguments are listed in their default order, which can easily be verified by using the standard way to get help about a function (in this case ``?mxMatrix()`` ).

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        myModel5     <- mxModel("model5",
                        mxMatrix(type="Full", nrow=1, ncol=1, values=4, name="Amatrix"),
                        mxAlgebra(expression=Amatrix+1, name="Bmatrix")
                        )
        myModel5run  <- mxRun(myModel5)

Note that all arguments are separated by commas.  In this case, we've also separated the arguments on different lines, but that is only for clarity.  No comma is needed after the last argument!  If you accidentally put one in, you get the generic error message *'argument is missing, with no default'* meaning that you forgot something and R doesn't know what it should be. The bracket on the following line closes the ``mxModel()`` statement.


Data functions
--------------

Most models will be fitted to data, not just a single number.  We will briefly introduce how to read data that are pre-packaged with the OpenMx library as well as reading in your own data.  All standard R utilities can be used here.  The critical part is to run an OpenMx model on these data, thus another OpenMx function ``mxData()`` is needed.

Reading Data
^^^^^^^^^^^^

The ``data`` function can be used to read sample data that has been pre-packaged into the OpenMx library. One such sample data set is called "demoOneFactor".  

.. cssclass:: input
..

   R Code
   
   .. code-block:: r

        data(demoOneFactor)

In order to read your own data, you will most likely use the ``read.table``, ``read.csv``, ``read.delim`` functions, or other specialized functions available from CRAN to read from 3rd party sources.  We recommend you install the package **psych** which provides succinct descriptive statistics with the ``describe()`` function.

.. cssclass:: input
..

   R Code
   
   .. code-block:: r
   
        require(psych)
        describe(demoOneFactor)

The output of this function is shown below.

    .. code-block:: r

           var   n  mean   sd median trimmed  mad   min  max range  skew kurtosis   se
        x1   1 500 -0.04 0.45  -0.03   -0.04 0.46 -1.54 1.22  2.77 -0.05     0.01 0.02
        x2   2 500 -0.05 0.54  -0.03   -0.04 0.55 -2.17 1.72  3.89 -0.14     0.05 0.02
        x3   3 500 -0.06 0.61  -0.03   -0.05 0.58 -2.29 1.83  4.12 -0.17     0.23 0.03
        x4   4 500 -0.06 0.73  -0.08   -0.05 0.75 -2.48 2.45  4.93 -0.08     0.05 0.03
        x5   5 500 -0.08 0.82  -0.08   -0.07 0.89 -2.62 2.18  4.80 -0.10    -0.23 0.04

Now that the data are accessible in R, we need to make them readable into our OpenMx model.

Data Source 
^^^^^^^^^^^

A ``mxData()`` function is used to construct a data source for the model.  OpenMx can handle fitting models to summary statistics and to raw data.

The most commonly used **summary statistics** are covariance matrices, means and correlation matrices; information on the variances is lost/unavailable with correlation matrices, so these are usually not recommended.

These days, the standard approach for model fitting applications is to use **raw data**, which is simply a data table or rectangular file with columns representing variables and rows representing subjects.  The primary benefit of this approach is that it handles datasets with missing values very conveniently and appropriately.

Covariance Matrix
^^^^^^^^^^^^^^^^^

We will start with an example using summary data, so we are specifying a covariance matrix by using the R function ``cov`` to generate a covariance matrix from the data frame.  In addition to reading in the actual covariance matrix as the first (``observed``) argument, we specify the ``type`` (one of "cov","cor","sscp" and "raw") and the number of observations (``numObs``).

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        exampleDataCov <- mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
    
We can view what *exampleDataCov* looks like for OpenMx.

    .. code-block:: r

         > 	exampleDataCov
         MxData 'data' 
         type : 'cov' 
         numObs : '500' 
         Data Frame or Matrix : 
                   x1        x2        x3        x4        x5
         x1 0.1985443 0.1999953 0.2311884 0.2783865 0.3155943
         x2 0.1999953 0.2916950 0.2924566 0.3515298 0.4019234
         x3 0.2311884 0.2924566 0.3740354 0.4061291 0.4573587
         x4 0.2783865 0.3515298 0.4061291 0.5332788 0.5610769
         x5 0.3155943 0.4019234 0.4573587 0.5610769 0.6703023
         Means : NA 
         Acov : NA 
         Thresholds : NA
    
Some models may include predictions for the mean(s).  We could add an additional ``means`` argument to the ``mxData`` statement to read in the means as well.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        exampleDataCovMeans <- mxData(observed=cov(demoOneFactor), 
                                   means=(colMeans(demoOneFactor), type="cov", numObs=500)
    
The output for *exampleDataCovMeans* would have the following extra lines.

    .. code-block:: r

        ....
        Means : 
                      x1          x2          x3          x4          x5
        [1,] -0.04007841 -0.04583873 -0.05588236 -0.05581416 -0.07555022
    
Raw Data
^^^^^^^^

Note that for most real life examples, raw data are the preferred option, except in cases where complete data are available on all variables included in the analyses.  In that situation, using summary statistics is faster.  To change the current example to use raw data, we would read in the data explicitly and specify the ``type`` as "raw".  The ``numObs`` is no longer required as the sample size is counted automatically.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        exampleDataRaw <- mxData(observed=demoOneFactor, type="raw")

Printing this MxData object would result in listing the whole data set.  We show just the first few lines here:

    .. code-block:: r

         > exampleData
         MxData 'data' 
         type : 'raw' 
         numObs : '500' 
         Data Frame or Matrix : 
                        x1            x2           x3           x4           x5
         1   -1.086832e-01 -0.4669377298 -0.177839881 -0.080931127 -0.070650263
         2   -1.464765e-01 -0.2782619339 -0.273882553 -0.154120074  0.092717293
         3   -6.399140e-01 -0.9295294042 -1.407963429 -1.588974090 -1.993461644
         4    2.150340e-02 -0.2552252972  0.097330513 -0.117444884 -0.380906486
         5    ....

The data to be used for our example are now ready in either **covariance matrix** or **raw data** format.

Model functions
---------------

We introduce here several new features by building a basic factor model to real data.  A useful tool to represent such a model is drawing a path diagram which is mathematically equivalent to equations describing the model.  If you're not familiar with the method of path analysis, we suggest you read one of the key reference books [LI1986]_.

.. [LI1986]  Li, C.C. (1986). Path Analysis - A Primer.  The Boxwood Press, Pacific Grove, CA.

Briefly, squares are used for observed variables; latent variables are drawn in circles.  One-headed arrows are drawn to represent causal relationships.  Correlations between variables are represented with two-headed arrows.  Double-headed paths are also used for variances of variables.  Below is a figure of a one factor model with five indicators (x1..x5). We have added a value of 1.0 to the variance of the latent variable **G** as a fixed value.  All the other paths in the models are considered free parameters and are to be estimated.

.. image:: graph/OneFactor.png
    :height: 2in
    
Variables
^^^^^^^^^

To specify this path diagram in OpenMx, we need to indicate which variables are observed or manifest and which are latent.  The ``mxModel()`` arguments ``manifestVars`` and ``latentVars`` both take a vector of variable names.   In this case the manifest variables are "x1", "x2", "x3", "x4", "x5" and the latent variable is "G".  The R function ``c()`` is used to build the vectors.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        manifests <- c("x1","x2","x3","x4","x5")
        latents <- c("G")
        
        manifestVars = manifests
        latentVars = latents

This could be written more succinctly as follows.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r
   
        manifestVars = names(demoOneFactor)
        latentVars = c("G")

because the R ``names()`` function call returns the vector of names that we want (the observed variables in the data frame "demoOneFactor").

Path Creation
^^^^^^^^^^^^^

Paths are created using the ``mxPath()`` function. Multiple paths can be created with a single invocation of the ``mxPath()`` function. 

- The ``from`` argument specifies the path sources, and the ``to`` argument specifies the path sinks.  If the ``to`` argument is missing, then it is assumed to be identical to the ``from`` argument. 
- The ``connect`` argument specifies the type of the source to sink connection, which can be one of five types.  For our example, we use the default "single" type in which the :math:`i^{th}` element of the ``from`` argument is matched with the :math:`i^{th}` element of the ``to`` argument, in order to create a path.  
- The ``arrows`` argument specifies whether the path is unidirectional (single-headed arrow, "1") or bidirectional (double-headed arrow, "2").  
- The next three arguments are vectors: ``free``, is a boolean vector that specifies whether a path is free or fixed; ``values`` is a numeric vector that specifies the starting value of the path; ``labels`` is a character vector that assigns a label to each free or fixed parameter.  Paths with the same labels are constrained to be equal, and OpenMx insists that paths equated in this way have the same fixed or free status; if this is not the case it will report an error.

To specify the path model above, we need to specify three different sets of paths.  The first are the single-headed arrows from the latent to the manifest variables, which we will put into the R object *causalPaths* as they represent causal paths.  The second set are the residuals on the manifest variables, referred to as *residualVars*.  The third ``mxPath()`` statement fixes the variance of the latent variable to one, and is called *factorVars*.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        causalPaths  <- mxPath(from=latents, to=manifests)
        residualVars <- mxPath(from=manifests, arrows=2)
        factorVars   <- mxPath(from=latents, arrows=2, free=FALSE, values=1.0)

Note that several arguments are optional.  For example, we omitted the ``free`` argument for *causalPaths* and *residualVars* because the default is 'TRUE' which applies in our example.  We also omitted the ``connect`` argument for all three paths.  The default "single" type automatically generates paths from every variable back to itself for all the variances, both the *residualVars* or the *factorVars*, as neither of those statements includes the ``to`` argument.  For the *causalPaths*, the default ``connect`` type will generate separate paths from the latent to each of the manifest variables.  To keep things simple, we did not include ``values`` or ``labels`` arguments as they are not strictly needed for this example, but this may not be true in general.  Once the variables and paths have been specified, the predicted covariance matrix will be generated from the implied path diagram in the backend of OpenMx using the RAM notation (see below).

Equations
^^^^^^^^^

For those more in tune with equations and matrix algebra, we can represent the model using matrix algebra rather than path specifications.  For reasons that may become clear later, the expression for the expected covariances between the manifest variables is given by  

.. math::
   :nowrap:

   \begin{eqnarray*} 
   \mbox{Cov} ( x_{ij}) = facLoadings * facVariances * facLoadings^\prime + resVariances
   \end{eqnarray*}

where *facLoadings* is a column vector of factor loadings, *facVariances* is a symmetric matrix of factor variances and *resVariances* is a diagonal matrix of residual variances.  You might have noticed the correspondence between *causalPaths* and *facLoadings*, between *residualVars* and *resVariances*, and between *factorVars* and *facVariances*.  To translate this model into OpenMx using the matrix specification, we will define the three matrices first using the ``mxMatrix()`` function, and then specify the algebra using the ``mxAlgebra()`` function.

Matrix Creation
^^^^^^^^^^^^^^^

The next three lines create three ``MxMatrix()`` objects, using the ``mxMatrix()`` function.  The first argument declares the ``type`` of the matrix, the second argument declares the number of rows in the matrix (``nrow``), and the third argument declares the number of columns (``ncol``).  The ``free`` argument specifies whether an element is a free or fixed parameter.  The ``values`` argument specifies the starting values for the elements in the matrix, and the ``name`` argument specifies the name of the matrix. 

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        mxFacLoadings  <-  mxMatrix(type="Full", nrow=5, ncol=1, 
                                    free=TRUE, values=0.2, name="facLoadings")
        mxFacVariances <-  mxMatrix(type="Symm", nrow=1, ncol=1, 
                                    free=FALSE, values=1, name="facVariances")
        mxResVariances <-  mxMatrix(type="Diag", nrow=5, ncol=5, 
                                    free=TRUE, values=1, name="resVariances")

Each ``MxMatrix()`` object is a container that stores five matrices of equal dimensions.  The five matrices stored in a ``MxMatrix()`` object are: ``free``, ``values``, ``labels``, ``lbound``, and ``ubound``.  ``Free`` stores a boolean vector that determines whether a element is free or fixed.  ``Values`` stores the current values of each element in the matrix.  ``Labels`` stores a character label for each element in the matrix. And ``lbound`` and ``ubound`` store the lower and upper bounds, respectively, for each element that is a free parameter.  If a element has no label, lower bound, or upper bound, then an NA value is stored in the element of the respective matrix.
 
Algebra Creation
^^^^^^^^^^^^^^^^

An ``mxAlgebra()`` function is used to construct an expression for any algebra, in this case the expected covariance algebra.  The first argument (``expression``) is the algebra expression that will be evaluated by the numerical optimizer.  The matrix operations and functions that are permitted in an MxAlgebra expression are listed in the help for the ``mxAlgebra`` function (obtained by ``?mxAlgebra``).  The algebra expression refers to entities according to ``name`` argument of the MxMatrix objects.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        mxExpCov     <- mxAlgebra(expression=facLoadings %*% facVariances %*% t(facLoadings) 
                                  + resVariances, name="expCov")

You can see a direct correspondence between the formula above and the expression used to create the expected covariance matrix *myExpCov*.

Expectation - Fit Function
--------------------------

To fit a model to data, the differences between the observed covariance matrix (the data, in this case the summary statistics) and model-implied expected covariance matrix are minimized using a fit function.  Fit functions are functions for which free parameter values are chosen such that the value of the fit function is minimized.  Now that we have specified data objects and path or matrix/algebra objects for the predicted covariances of our model, we need to link the two and execute them which is typically done with ``mxExpectation()`` and ``mxFitFunction()`` statements.  PS. These two statements replace the ``mxObjective()`` functions`` in earlier versions of OpenMx.  

RAM Expectation 
^^^^^^^^^^^^^^^

When using a path specification of the model, the fit function is always ``RAM`` which is indicated by using the ``type`` argument.  We don't have to specify the fit function explicitly with an ``mxExpectation()`` and ``FitFunction()`` argument, instead we simply add the following argument to the model.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        type="RAM"
    
To gain a better understanding of the RAM principles, we recommend reading [RAM1990]_

Normal Expectation
^^^^^^^^^^^^^^^^^^

When using a matrix specification, ``mxExpectationNormal()`` defines how model expectations are calculated using the matrices/algebra implied by the ``covariance`` argument and optionally the ``means``.  For this example, we are specifying an expected covariance algebra (``covariance``) omitting an expected means algebra.  The expected covariance algebra is referenced according to its name, i.e. the ``name`` argument of the MxAlgebra created above.  We also need to assign ``dimnames`` for the rows and columns of this covariance matrix, such that a correspondence can be determined between the expected and the observed covariance matrices.  Subsequently we are specifying a maximum likelihood fit function with the ``mxFitFunctionML()`` statement.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        expectCov    <- mxExpectationNormal(covariance="expCov", 
                                            dimnames=names(demoOneFactor))
        funML        <- mxFitFunctionML()

The above expectation and fit function can be used when fitting to covariance matrices.  A model for the predicted means is optional.  However, when fitting to raw data, an expectation has to be used that specifies both a model for the means and for the covariance matrices, paired with the appropriate fit function.  In the case of raw data, the ``mxFitFunctionML()`` function uses full-information maximum likelihood to provide maximum likelihood estimates of free parameters in the algebra defined by the ``covariance`` and ``means`` arguments. The ``covariance`` argument takes an ``MxMatrix`` or ``MxAlgebra`` object, which defines the expected covariance of an associated ``MxData`` object. Similarly, the ``means`` argument takes an ``MxMatrix`` or ``MxAlgebra`` object to define the expected means of an associated ``MxData`` object. The ``dimnames`` arguments takes an optional character vector. This vector is assigned to be the ``dimnames`` of the means vector, and the row and columns ``dimnames`` of the covariance matrix. 

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        expectCovMeans <- mxExpectationNormal(covariance="expCov", means="expMeans", 
                                              dimnames=names(demoOneFactor))
        funML        <- mxFitFunctionML()

Raw data can come in two forms, continuous or categorical.  While **continuous data** have an unlimited number of possible values, their frequencies typically form a normal distribution.

There are basically two flavors of **categorical data**.  If only two response categories exist, for example Yes and No, or affected and unaffected, we are dealing with binary data.  Variables with three or more ordered categories are considered ordinal.

Continuous Data
^^^^^^^^^^^^^^^

When the data to be analyzed are continuous, and models are fitted to raw data, the ``mxFitFunctionML()`` function will take two arguments, the ``covariance`` and the ``means`` argument, as well as ``dimnames`` to match them up with the observed data.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        expectRaw    <- mxExpectationNormal(covariance="expCov", means="expMeans", 
                                            dimnames=manifests)
        funML        <- mxFitFunctionML()

If the variables to be analyzed have at least 15 possible values, we recommend to treat them as continuous data.  As will be discussed later in the documentation, the power of the study is typically higher when dealing with continuous rather than categorical data.

Categorical Data
^^^^^^^^^^^^^^^^

For categorical - be they binary or ordinal - data, an additional argument is needed for the ``mxFitFunctionML()`` function, besides the ``covariance`` and ``means`` arguments, namely the ``thresholds`` argument.
    
.. cssclass:: input
..

    OpenMx Code

    .. code-block:: r

        expFunOrd    <- mxExpectationNormal(covariance="expCov", means="expMeans", 
                                            thresholds="expThres", dimnames=manifests)
        funML        <- mxFitFunctionML()

For now, we will stick with the factor model example and fit it to covariance matrices, calculated from the raw continuous data.


Methods
-------

We have introduced two ways to create a model.  One is the **path method**, in which observed and latent variables are specified as well as the causal and correlational paths that connect the variables to form the model.  This method may be more intuitive as the model maps on directly to the diagram.  This of course assumes that the path diagram is drawn mathematically correct.  Once the model is 'drawn' or specified correctly in this way, OpenMx translates the paths into RAM notation for the predicted covariance matrices.

Alternatively, we can specify the model using the **matrix method** by creating the necessary matrices and combining them using algebra to generate the expected covariance matrices (and optionally the mean/threshold vectors).  Although less intuitive, this method provides greater flexibility for developing more complex models.  Let us look at examples of both.

Path Method
^^^^^^^^^^^

We have previously generated all the pieces that go into the model, using the path method specification.  As we have discussed before, the ``mxModel()`` function is somewhat of a swiss-army knife.  The first argument to the ``mxModel()`` function can be an argument of type ``name`` (and appear in quotes), in which case it is a newly generated model, or it can be a previously defined model object.  In the latter case, the new model 'inherits' all the characteristics (arguments) of the old model, which can be changed with additional arguments.  An ``mxModel()`` can contain ``mxData()``, ``mxPath()``, ``mxExpectation()``, ``mxFitFunction`` and other ``mxModel()`` statements as arguments.

The  following ``mxModel()`` function is used to create the 'one-factor' model, shown on the path diagram above.  The first argument is a ``name``, thus we are specifying a new model, called "One Factor".  By specifying the ``type`` argument to equal "RAM", we create a path style model. A RAM style model must include a vector of manifest variables (``manifestVars=``) and a vector of latent variables (``latentVars=``).   We then include the arguments for reading the example data *exampleDataCov*, and those that specify the paths of the path model *causalPaths*, *residualVars*, and *factorVars* which we created previously.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        factorModel1 <- mxModel(name="One Factor", 
            type="RAM",
            manifestVars=manifests,
            latentVars=latents,
            exampleDataCov, causalPaths, residualVars, factorVars)

When we display the contents of this model, note that we now have manifest and latent variables specified.  By using ``type``="RAM" we automatically use the expectation ``mxExpectationRAM`` which translates the path model into RAM specification [RAM1990] as reflected in the matrices **A**, **S** and **F**,  and the function ``mxFitFunctionML()``.  Briefly, the **A** matrix contains the asymmetric paths, which are the unidirectional paths in the *causalPaths* object, and represent the factor loadings from the latent variable onto the manifest variables.  The **S** matrix contains the symmetric paths which include both the bidirectional paths in *residualVars* and in *factorVars*.  The **F** matrix is the filter matrix.

The formula :math:`F(I-A)^-1*S*(I-A)^-1'F'`, where I is an identity matrix, :math:`^-1` denotes the inverse and ' the transpose, generates the expected covariance matrix.

   .. code-block:: r

        >     factorModel1
        MxModel 'One Factor' 
        type : RAM 
        $matrices : 'A', 'S', and 'F' 
        $algebras :  
        $constraints :  
        $intervals :  
        $latentVars : 'G' 
        $manifestVars : 'x1', 'x2', 'x3', 'x4', and 'x5' 
        $data : 5 x 5 
        $data means : NA
        $data type: 'cov' 
        $submodels :  
        $expectation : MxExpectationRAM 
        $fitfunction : MxFitFunctionML 
        $compute : NULL 
        $independent : FALSE 
        $options :  
        $output : FALSE 

You can verify that after running the model, the new R object *factorFit* has similar arguments, except that they now contain the estimates from the model rather than the starting values.  For example, we can look at the values in the **A** matrix in the built model *factorModel*, and in the fitted model *factorFit*.  We will get back to this later.  Note also that from here on out, we use the convention the R object containing the built model will end with *Model* while the R object containing the fitted model will end with *Fit*.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        factorFit1 <- mxRun(factorModel1)

We can inspect the values of the **A** matrix in *factorModel1* and *factorFit1* respectively as follows.

    .. code-block:: r

        > factorModel1$A$values
           x1 x2 x3 x4 x5 G
        x1  0  0  0  0  0 0
        x2  0  0  0  0  0 0
        x3  0  0  0  0  0 0
        x4  0  0  0  0  0 0
        x5  0  0  0  0  0 0
        G   0  0  0  0  0 0

        > factorFit1$A$values 
           x1 x2 x3 x4 x5         G
        x1  0  0  0  0  0 0.3971521
        x2  0  0  0  0  0 0.5036611
        x3  0  0  0  0  0 0.5772414
        x4  0  0  0  0  0 0.7027737
        x5  0  0  0  0  0 0.7962500
        G   0  0  0  0  0 0.0000000

We can also specify all the arguments directly within the ``mxModel()`` function, using the **classical** style, as follows.  The script reads data from disk, creates the one factor model, fits the model to the observed covariances, and prints a summary of the results. 

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        data(demoOneFactor)
        manifests <- names(demoOneFactor)
        latents   <- c("G")
        
        factorModel1 <- mxModel(name="One Factor", 
            type="RAM",
            manifestVars=manifests,
            latentVars=latents,
            mxPath(from=latents, to=manifests),
            mxPath(from=manifests, arrows=2),
            mxPath(from=latents, arrows=2, free=FALSE, values=1.0), 
            mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
        )
        
        factorFit1 <- mxRun(factorModel1)
        summary(factorFit1)
    
For more details about the summary and alternative options to display model results, see below.

Matrix Method
^^^^^^^^^^^^^

We will now re-create the model from the previous section, but this time we will use a matrix specification technique. The script reads data from disk, creates the one factor model, fits the model to the observed covariances, and prints a summary of the results. 

We have already created separate objects for each of the parts of the model, which can then be combined in an ``mxModel()`` statement at the end.  To repeat ourselves, the name of an OpenMx entity bears no relation to the R object that is used to identify the entity. In our example, the object "mxFacLoadings" stores a value that is a MxMatrix object with the name "facLoadings".  Note, however, that it is not necessary to use different names for the name within the ``mxMatrix`` object and the name of the R object generated with the statement.  For more complicated models, using the same name for both rather different entities, may make it easier to keep track of the various pieces.  For now, we will use different names to highlight which one should be used in which context.
 
.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        data(demoOneFactor)
        
        factorModel2 <- mxModel(name="One Factor",
            exampleDataCov, mxFacLoadings, mxFacVariances, mxResVariances, 
            mxExpCov, expectCov, funML)
        factorFit2 <- mxRun(factorModel2)
        summary(factorFit2)

Alternatively, we can write the script in the **classical** style and specify  all the matrices, algebras, objective function and data as arguments to the ``mxModel()``.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        data(demoOneFactor)
        
        factorModel2 <- mxModel(name="One Factor",
            mxMatrix(type="Full", nrow=5, ncol=1, free=TRUE, values=0.2, name="facLoadings"),
            mxMatrix(type="Symm", nrow=1, ncol=1, free=FALSE, values=1, name="facVariances"),
            mxMatrix(type="Diag", nrow=5, ncol=5, free=TRUE, values=1, name="resVariances"),
            mxAlgebra(expression=facLoadings %*% facVariances %*% t(facLoadings) 
                                + resVariances, name="expCov"),
            mxExpectationNormal(covariance="expCov", dimnames=names(demoOneFactor)),
            mxFitFunctionML()
            mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
        )
        
        factorFit2 <- mxRun(factorModel2)
        summary(factorFit2)

Now that we've specified the model with both methods, we can run both examples and verify that they indeed provide the same answer by inspecting the two fitted R objects *factorFit1* and *factorFit2*.

Output
------

We can generate output in a variety of ways.  As you might expect, the **summary** function summarizes the model, including data, model parameters, goodness-of-fit and run statistics.

Note that the fitted model is an R object that can be further manipulated, for example, to output specific parts of the model or to use it as a basis for developing an alternative model.

Model Summary
^^^^^^^^^^^^^

The summary function (``summary(modelname)``) is a convenient method for displaying the highlights of a model after it has been executed.  Many R functions have an associated ``summary()`` function which summarizes all key aspects of the model.  In the case of OpenMx, the ``summary(model)`` includes a summary of the data, a list of all the free parameters with their name, matrix element locators, parameter estimate and standard error, as well as lower and upper bounds if those were assigned.  Currently the list of goodness-of-fit statistics printed include the number of observed statistics, the number of estimated parameters, the degrees of freedom, minus twice the log-likelihood of the data, the number of observations, the chi-square and associated p-value and several information criteria.  Various time-stamps and the OpenMx version number are also displayed.

   .. code-block:: r

        >     summary(factorFit1)
        data:
        $`One Factor.data`
        $`One Factor.data`$cov
                  x1        x2        x3        x4        x5
        x1 0.1985443 0.1999953 0.2311884 0.2783865 0.3155943
        x2 0.1999953 0.2916950 0.2924566 0.3515298 0.4019234
        x3 0.2311884 0.2924566 0.3740354 0.4061291 0.4573587
        x4 0.2783865 0.3515298 0.4061291 0.5332788 0.5610769
        x5 0.3155943 0.4019234 0.4573587 0.5610769 0.6703023


        free parameters:
          name matrix row col   Estimate   Std.Error Std.Estimate      Std.SE lbound ubound
        1  One Factor.A[1,6]      A  x1   G 0.39715182 0.015549708   0.89130932 0.034897484              
        2  One Factor.A[2,6]      A  x2   G 0.50366066 0.018232433   0.93255458 0.033758321              
        3  One Factor.A[3,6]      A  x3   G 0.57724092 0.020448313   0.94384664 0.033435037              
        4  One Factor.A[4,6]      A  x4   G 0.70277323 0.024011318   0.96236250 0.032880581              
        5  One Factor.A[5,6]      A  x5   G 0.79624935 0.026669339   0.97255562 0.032574489              
        6  One Factor.S[1,1]      S  x1  x1 0.04081418 0.002812716   0.20556770 0.014166734              
        7  One Factor.S[2,2]      S  x2  x2 0.03801997 0.002805791   0.13034196 0.009618951              
        8  One Factor.S[3,3]      S  x3  x3 0.04082716 0.003152305   0.10915353 0.008427851              
        9  One Factor.S[4,4]      S  x4  x4 0.03938701 0.003408870   0.07385841 0.006392303              
        10 One Factor.S[5,5]      S  x5  x5 0.03628708 0.003678556   0.05413557 0.005487924              

        observed statistics:  15 
        estimated parameters:  10 
        degrees of freedom:  5 
        -2 log likelihood:  -3648.281 
        saturated -2 log likelihood:  -3655.665 
        number of observations:  500 
        chi-square:  7.384002 
        p:  0.1936117 
        Information Criteria: 
             df Penalty Parameters Penalty Sample-Size Adjusted
        AIC:  -2.615998           27.38400                   NA
        BIC: -23.689038           69.53008             37.78947
        CFI: 0.9993583 
        TLI: 0.9987166 
        RMSEA:  0.03088043 
        timestamp: 2014-04-10 10:23:07 
        frontend time: 0.02934313 secs 
        backend time: 0.005492926 secs 
        independent submodels time: 1.907349e-05 secs 
        wall clock time: 0.03485513 secs 
        cpu time: 0.03485513 secs 
        openmx version number: 999.0.0

The table of free parameters requires a little more explanation.  First, ``<NA>`` is given for the name of elements that were not assigned a label.  Second, the columns 'row' and 'col' display the variables at the tail of the paths and the variables at the head of the paths respectively.  Third, standard errors are calculated.  We will discuss the use of standard errors versus confidence intervals later on.

Model Evaluation 
^^^^^^^^^^^^^^^^

The ``mxEval()`` function should be your primary tool for observing and manipulating the final values stored within a MxModel object.  The simplest form of the ``mxEval()`` function takes two arguments: an ``expression`` and a ``model``. The expression can be **any** arbitrary expresssion to be evaluated in R.  That expression is evaluated, but the catch is that any named entities or parameter names are replaced with their current values from the model.  The model can be either a built or a fitted model.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        myModel6      <- mxModel('topmodel', 
            mxMatrix('Full', 1, 1, values=1, free=TRUE, labels='p1', name='A'),
            mxModel('submodel', 
                mxMatrix('Full', 1, 1, values=2, free=FALSE, labels='p2', name='B')
            )
        )
        myModel6Run   <- mxRun(myModel6)

The example above has a model ("submodel") embedded in another model ("topmodel").  Note that the name of the arguments can be omitted if they are used in the default order (``type``, ``nrow`` and ``ncol``).

The ``expression`` of the ``mxEval`` statement can include both matrices, algebras as well as matrix element labels, each taking on the value of the model specified in the ``model`` argument.  To reinforce an earlier point, it is not necessary to restrict the expression only to valid MxAlgebra expressions.  In the following example, we use the ``harmonic.mean()`` function from the ``psych`` package.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r
 
        mxEval(A + submodel.B + p1 + p2, myModel6)       # initial values
        mxEval(A + submodel.B + p1 + p2, myModel6Run)    # final values

        library(psych)
        nVars <- 4
        mxEval(nVars * harmonic.mean(c(A, submodel.B)), myModel6)

When the name of an entity in a model collides with the name of a built-in or user-defined function in R, the named entity will supercede the function.  We strongly advice against naming entities with the same name as the predefined functions or values in R, such as `c`, `T`, and `F` among others.

The ``mxEval()`` function allows the user to inspect the values of named entities without explicitly poking at the internals of the components of a model.  We encourage the use of ``mxEval()`` to look at the state of a model either before the execution of a model or after model execution.


Indexing Operator
^^^^^^^^^^^^^^^^^

MxModel objects support the ``$`` operator, also known as the list indexing operator, to access all the components contained within a model.  Here is an example collection of models that will help explain the uses of the ``$`` operator:

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r
   
        myModel7 <- 
            mxModel('topmodel', 
                mxMatrix(type='Full', nrow=1, ncol=1, name='A'),
                mxAlgebra(A, name='B'),
                mxModel('submodel1', 
                    mxConstraint(topmodel1.A == topmodel1.B, name = 'C'),
                    mxModel('undersub1', mxData(diag(3), type='cov', numObs=10) )
                ),
                mxModel('submodel2', 
                    mxFitFunctionAlgebra('topmodel1.A')
                )
            )

Access Elements
^^^^^^^^^^^^^^^

The first useful trick is entering the string ``model$`` in the R interpreter and then pressing the TAB key.  You should see a list of all the named entities contained within the ``model`` object.

    .. code-block:: r

        > model$
        model$A                    
        model$B                    
        model$submodel1
        model$submodel2            
        model$submodel1.C          
        model$undersub1
        model$undersub1.data
        model$submodel2.fitfunction

The named entities of the model are displayed in one of three modes. 

#. All of the submodels contained within the parent model are accessed by using their unique model name (``submodel1``, ``submodel2``, and ``undersub1``).  

#. All of the named entities contained within the parent model are displayed by their names (``A`` and ``B``).  

#. All of the named entities contained by the submodels are displayed in the ``modelname.entityname`` format (``submodel1.C``, ``submodel2.objective``, and ``undersub1.data``). 

Modify Elements
^^^^^^^^^^^^^^^

The list indexing operator can also be used to modify the components of an existing model. There are three modes of using the list indexing operator to perform modifications, and they correspond to the three modes for accessing elements.

In the first mode, a submodel can be replaced using the unique name of the submodel or even eliminated.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r

        # replace 'submodel1' with the contents of the mxModel() expression
        model$submodel1 <- mxModel(...)      
        # eliminate 'undersub1' and all children models
        model$undersub1 <- NULL              

In the second mode, the named entities of the parent model are modified using their names.  Existing matrices can be eliminated or new matrices can be created.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r
   
        # eliminate matrix 'A'
        model$A <- NULL
        # create matrix 'D'
        model$D <- mxMatrix(...)             

In the third mode, named entities of a submodel can be modified using the ``modelname.entityname`` format.  Again existing elements can be eliminated or new elements can be created.

.. cssclass:: input
..

   OpenMx Code
   
   .. code-block:: r
   
        # eliminate constraint 'C' from submodel1
        model$submodel1.C <- NULL
        # create algebra 'D' in undersub1
        model$undersub1.D <- mxAlgebra(...)         
        # create 'undersub2' as a child model of submodel1
        model$submodel1.undersub2 <- mxModel(...)   

Keep in mind that when using the list indexing operator to modify a named entity within a model, the name of the created or modified entity is always the name on the left-hand side of the ``<-`` operator.  This feature can be convenient, as it avoids the need to specify a name of the entity on the right-hand side of the ``<-`` operator.


Classes
-------

We have introduced a number of OpenMx functions which correspond to specific classes which are summarized below. 
The basic unit of abstraction in the OpenMx library is the model.  A model serves as a container for a collection of matrices, algebras, constraints, expectation, fit functions, data sources, and nested sub-models.  In the parlance of R, a model is a value that belongs to the class MxModel that has been defined by the OpenMx library.  The following table indicates what classes are defined by the OpenMx library.

                    +--------------------+---------------------+
                    | entity             | S4 class            |
                    +====================+=====================+
                    | model              | MxModel             | 
                    +--------------------+---------------------+
                    | data source        | MxData              |
                    +--------------------+---------------------+
                    | matrix             | MxMatrix            |
                    +--------------------+---------------------+
                    | algebra            | MxAlgebra           |
                    +--------------------+---------------------+
                    | expectation        | MxExpectationRAM    |
                    |                    | MxExpectationNormal |
                    +--------------------+---------------------+
                    | fit function       | MxFitFunctionML     |
                    +--------------------+---------------------+                    
                    | constraint         | MxConstraint        |
                    +--------------------+---------------------+

All of the entities listed in the table are identified by the OpenMx library by the name assigned to them.  A name is any character string that does not contain the "." character.  In the parlance of the OpenMx library, a model is a container of named entities.  The name of an OpenMx entity bears no relation to the R object that is used to identify the entity. In our example, the object ``factorModel`` is created with the ``mxModel()`` function and stores a value that is an "MxModel" object with the name 'One Factor'.

.. [RAM1990]  McArdle, J.J. & Boker, S.M. (1990). RAMpath: Path diagram software. Denver: Data Transforms Inc.
