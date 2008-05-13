###########################################################################/**
# @RdocClass MxJob
#
# @title "The MxJob Class"
#
# \description{
#
#  The MxJob is the basic unit of computation.
# 
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{model}{A MxMatrix object.}
#   \item{objective}{A MxObjective object.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
#
#*/###########################################################################
setConstructorS3("MxJob", function(model, objective) {

   if (missing(model)) model <- NA;
   if (missing(objective)) objective <- NA;


   extend(Object(), "MxJob",
	.model=model,
	.objective=objective);

})

#########################################################################/**
# @RdocMethod createMxClosure
#
# @title "Create a MxJob Closure"
# 
# \description{
#    Create a MxJob Closure that will perform optimization.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxJob object.}
#  \item{use_R}{A boolean value whether or not to use R in computation.}
#  \item{...}{Unused.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("createMxClosure", "MxJob", function(this, use_R = FALSE, ...) {
   if (use_R) {
      createMxJobClosureR(this$.objective, this);
   } else {
      createMxJobClosureC(this$.objective, this);
   }
});

