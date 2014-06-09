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

MxModels are created in R, consisting of R data structures. An MxModel
can only have 1 MxData and 1 MxExpectation.  This implies that
multigroup and multilevel models require more than 1 MxModel.
MxRun translates all the information contained in an R
MxModel into corresponding C data structures. At this stage, as much
as possible is checked to make sure that the model is correctly
specified to avoid errors during optimization.

* MxMatrices are torn apart. Free variable labels are matched to
  create a list of free variables recording locations where they are
  stored in which matrices. Starting values and upper and lower bounds
  are collected for each free variable.
* The model tree is flattened. The hierarchy is preserved in
  MxBaseExpectation container and submodel slots, but otherwise, the
  whole model tree is flattened into lists of matrices, algebras,
  expectations, and fitfunctions. This treatment is applied separately
  to MxModels marked as independent.
* MxCompute controls what happens in the backend.
* After the backend returns, the model is updated with everything that
  changed.
* The model\@runstate slot is suppose to preserve the post-backend
  state of the model so that summary output doesn't change after the
  user changes something in the model after mxRun.

Note that MxEval is almost pure R code.

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
the granularity of machines or processes with a machine while OpenMP
is only concerned with exploiting all the CPU within a single machine.
