OpenMx Style Guide
==================

The goal of the OpenMx Style Guide is to make our OpenMx scripts easier to read, share, and verify.

The '$' Operator
----------------

MxModel objects support the '$' operator, also known as the list indexing operator, to access all the components contained within a model.  Here is an example collection of models that will help explain the uses of the '$' operator:

.. code-block:: r

    model <- mxModel('topmodel', 
        mxMatrix(type='Full', nrow=1, ncol=1, name='A'),
        mxAlgebra(A, name='B'),
        mxModel('submodel1', 
            mxConstraint(topmodel1.A == topmodel1.B, name = 'C'),
            mxModel('undersub1', 
                mxData(diag(3), type='cov', numObs=10)
            )
        ),
        mxModel('submodel2', 
            mxAlgebraObjective('topmodel1.A')
        )
    )

Accessing Elements
^^^^^^^^^^^^^^^^^^

The first useful trick is entering the string ``model$`` in the R interpreter and then pressing the TAB key.  You should see a list of all the named entities contained within the ``model`` object::

   > model$
   model$A                    model$submodel2
   model$B                    model$submodel2.objective
   model$submodel1            model$undersub1
   model$submodel1.C          model$undersub1.data

The named entities of the model are displayed in one of three modes. In the first mode, all of the submodels contained within the parent model are accessed by using their unique model name (``submodel1``, ``submodel2``, and ``undersub1``).  In the second mode, all of the named entities contained within the parent model are displayed by their names (``A`` and ``B``).  In the third mode, all of the named entities contained by the submodels are displayed in the ``modelname.entityname`` format (``submodel1.C``, ``submodel2.objective``, and ``undersub1.data``). The three models will become even more important in the next section, *Modifying Elements*, so make sure you are comfortable with them before moving on.

Modifying Elements
^^^^^^^^^^^^^^^^^^

The list indexing operator can also be used to modify the components of an existing model. There are three modes of using the list indexing operator to perform modifications, and they correspond to the three models for accessing elements.

In the first mode, a submodel can be replaced using the unique name of the submodel, or even eliminated::

    # eliminate 'undersub1' and all children models
    model$undersub1 <- NULL              
    # replace 'submodel1' with the contents of the mxModel() expression
    model$submodel1 <- mxModel(...)      

In the second mode, the named entities of the parent model are modified using their names::

    # eliminate matrix 'A'
    model$A <- NULL                      
    # create matrix 'D'
    model$D <- mxMatrix(...)             

In the third mode, named entities of a submodel can be modified using the ``modelname.entityname`` format::

    # eliminate constraint 'C' from submodel1
    model$submodel1.C <- NULL                   
    # create algebra 'D' in undersub1
    model$undersub1.D <- mxAlgebra(...)         
    # create 'undersub2' as a child model of submodel1
    model$submodel1.undersub2 <- mxModel(...)   

Keep in mind that when using the list indexing operator to modify a named entity within a model, the name of the created or modified entity is always the name on the left-hand side of the ``<-`` operator.  This feature can be convenient, as it avoids the need to specify a name of the entity on the right-hand side of the ``<-`` operator.

mxEval()
--------

The ``mxEval()`` function should be your primary tool for observing and manipulating the final values stored within a MxModel object.  The simplest form of the mxEval function takes two arguments: an expression and a model. The expression can be **any** arbitrary expresssion to be evaluated in R.  That expression is evaluated, but the catch is that any named entities or parameter names are replaced with their current values from the model.

.. code-block:: r

    model <- mxModel('topmodel', 
        mxMatrix('Full', 1, 1, values=1, free=TRUE, labels='p1', name='A'),
        mxModel('submodel1', 
            mxMatrix('Full', 1, 1, values=2, free=FALSE, labels='p2', name='B')
        ),
        mxModel('submodel2', 
            mxMatrix('Full', 1, 1, values=3, name = 'B')
        )
    )

   modelOut <- mxRun(model)
   mxEval(A + submodel1.B + submodel2.B + p1 + p2, model)       # initial values
   mxEval(A + submodel1.B + submodel2.B + p1 + p2, modelOut)    # final values

To reinforce an earlier point, it is not necessary to restrict the expression to only valid MxAlgebra expressions.  In the following example, we use the ``harmonic.mean`` function from the psych package.

.. code-block:: r

   library(psych)
   nVars <- 3
   mxEval(nVars * harmonic.mean(c(A, submodel1.B, submodel2.B)), model)

When the name of an entity in a model collides with the name of a built-in or user-defined function in R, the named entity will supercede the function.  We strongly advice against naming entities with the same name as the predefined functions or values in R, such as `c`, `T`, and `F` among others.

The ``mxEval()`` function allows the user to inspect the values of named entities without explicitly poking at the internals of the components of a model.  We encourage the use of mxEval() to look at the state of a model either before the execution of a model or after execution.
