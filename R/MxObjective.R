###########################################################################/**
# @RdocClass MxObjective
#
# @title "The MxObjective Class"
#
# \description{
#
#  This is the abstract superclass of all objective function types.
# 
#  @classhierarchy
# }
# 
# @synopsis
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
#
#*/###########################################################################
setConstructorS3("MxObjective", function() {

  extend(Object(), "MxObjective");

});

#########################################################################/**
# @RdocMethod createMxJobClosureR
#
# @title "Create a Closure in R"
# 
# \description{
#    Create a MxJob Closure that will perform optimization in R.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxObjective object.}
#  \item{job}{The MxJob object.}
#  \item{...}{Unused.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("createMxJobClosureR", "MxObjective", function(this, job, ...) {});



#########################################################################/**
# @RdocMethod createMxJobClosureC
#
# @title "Create a Closure in C"
# 
# \description{
#    Create a MxJob Closure that will perform optimization in C.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxObjective object.}
#  \item{job}{The MxJob object.}
#  \item{...}{Unused.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("createMxJobClosureC", "MxObjective", function(this, job, ...) {});