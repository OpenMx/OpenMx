###########################################################################/**
# @RdocClass CovarianceObjective
#
# @title "The CovarianceObjective Class"
#
# \description{
#
#  This is covariance objective function type.
# 
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#  \item{expected}{The expected covariance algebra (MxAlgebra object).}
#  \item{observed}{The observed covariance matrix.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
#
#*/###########################################################################
setConstructorS3("CovarianceObjective", function(expected, observed) {

   if (missing(expected)) expected <- NA;
   if (missing(observed)) observed <- NA;


   extend(MxObjective(), "CovarianceObjective",
	.expected=expected,
	.observed=observed);

})

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
setMethodS3("createMxJobClosureR", "CovarianceObjective", function(this, job, ...) {
      
   model <- job$.model;
   model$makeFreeParametersList(); 

   function() {

      observedCov <- this$.observed; # this is a matrix
      expectedCov <- this$.expected; # this is an MxAlgebra statement
      startValues <- model$getFreeParameters();
      expectedCov$translateAlgebra();

      objectiveFunction <- function(par) {
         model$updateMatrices(par);
         expectedCov$setDirtyBit(TRUE);
         predictedCov <- expectedCov$evaluateTranslation();
         functionValue <- sum(observedCov * solve(predictedCov)) + log(det(predictedCov));
         return(functionValue);
      }

   # Optimize
   outNLM <- nlm(objectiveFunction, 
   		startValues, 
   		hessian=TRUE, 
   		iterlim=1000, 
   		print.level=2,
   		typsize=abs(startValues));

   }

})
