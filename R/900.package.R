#########################################################################/**
# @RdocPackage OpenMx
#
# \description{
#     The OpenMx Project intends to rewrite and extend the popular
#     statistical package Mx to address the challenges facing a large
#     range of modern statistical problems such as:
# \itemize{
#    \item The difficulty of measuring behavioral traits
#    \item The availability of technologies - such as such as
#          magnetic resonance imaging, continuous physiological
#          monitoring and microarrays - which generate extremely
#          large amounts of data often with complex time-dependent
#          patterning, and
#    \item Increased sophistication in the statistical models used
#          to analyze the data.
# }
#    To address these problems, OpenMx will rewrite the Mx Structural
#    Equation Modeling software so as to be:
# \itemize{
#    \item Split into modules that interoperate with the R statistical package,
#    \item Released as open source so as to provide a stable path for
#           future maintenance and development, and
#    \item Integrated with the VDL parallel workflow software.
#  }
#    Grid/parallel computing and data management using VDL will provide
#    significant speedup for processing large (up to multi-terabyte)
#    data sets, through the use of analytical workflows that provide
#    detailed provenance tracking and annotation of derived results.
#    Revised algorithms for model fitting and optimization will increase
#    both the scope of the software and its performance. Both the code
#    and its use will be documented and disseminated at
#    national and international workshops.
# }
# \section{To get started}{
#   To get started, see:
#   \enumerate{
#     \item @see "MxModel" - MxModel is a collection of
#          MxMatrix and MxAlgebra objects and the associated operations
#          that can be performed on these objects.
#     \item @see "MxMatrix" - The abstract superclass of all MxMatrix types.
#     \item @see "MxAlgebra" - The class stores a matrix algebra
#     expression to be evaluated at some later point.
#     \item @see "MxObjective" - The is the abstract superclass of
#     all objective function types.
#     \item @see "MxJob" - The MxJob is the basic unit of computation.
#     \item \code{demo()} - Watch the demos for the OpenMx package.
#   }
# }
#
#
#*/#########################################################################
