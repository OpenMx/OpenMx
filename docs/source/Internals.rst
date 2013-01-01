OpenMx Internal Architecture
============================

OpenMx is a language plugin designed to perform linear algebra.  OpenMx
is too limited to be seen as a general purpose programming language,
but it is also more expressive than a typical C library.  OpenMx was
developed as an R package. It borrows R's syntax parser to encode
algebra expressions and the internals are rich with the R procedures
needed to smoothly integrate with R. However, it is not inconceivable
that OpenMx could eventually be split into a pure C backend with an R
or Python or Julia front end.

The fundamental data structure in OpenMx is the 2 dimensional matrix
of double precision real numbers. This is inflexible. Data of
dimensionality greater than 2 is cumbersome to model. However, the
benefit of this inflexibility is speed. Internal operators can make
assumptions and go faster because matrices all have the same
structure.

Algebras are reified functions of matrices. Information flows through
algebras, in one direction, during model fitting. Algebras are not
compiled to assembly. Algebra use a simple byte code language to
represent algebraic operators. Dependencies between algebras are
recorded such that recalculation is kept to a minimum.  What
transformations cannot easily be expressed in algebras can be put into
algebra-like objects such as MxExpectation or MxFitFunction.

MxModel Lifecycle
-----------------

MxEval will evaluate some piece of a model, however, the real action
happens in MxRun. MxModels are created in R, consisting of R data
structures. MxRun translates all the information contained in an R
MxModel into corresponding C data structures. At this stage, as much
as possible is checked to make sure that the model is correctly
specified to avoid errors during optimization. For MxExpectation and
MxFitFunction, the initFun is called. Algebras have a chance to check
conformability.  A list of free variables is created recording
locations where they are stored in which matrices. The optimizer is
invoked. Free variables are moved around, repopulateFun is invoked,
and everything is recomputed according to free variable
dependencies. MxFitFunction.gradient can speed this up.  Eventually,
the optimizer will become satisfied (or give up). Thereafter, the free
variables will be further probed to compute the Hessian, confidence
intervals, and standard errors (also see
getStandardErrorFun). populateAttrFun and setFinalReturns will see to
it that all the freshly baked numbers are copied back into R data
structures. All C-side memory is freed.

Parallelization Facilities
--------------------------

There 2 facilities for OpenMx to take advantage of parallelism. The
first facility is OpenMP. At the time of writing, only the
FitFunctions FIMLFitFunction and RowFitFunction use OpenMP. To allow
threads to work without interfering with each other, both of these
MxFitFunctions request that the complete model (C data structures) be
copied such that each thread has its own copy. This model copying is
optional. With adequate care, it is also possible to write OpenMP
optimizations without the added complexity of duplicating the model
(see sadmvn.f).  The second facility is Snowfall. Snowfall works at
the granularity of machines while OpenMP is only concerned with
exploiting all the CPU within a single machine.

Questions
---------

What is the flat model?

How do submodels work?
