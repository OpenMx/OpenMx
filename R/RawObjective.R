###########################################################################/**
# @RdocClass RawObjective
#
# @title "The RawObjective Class"
#
# \description{
#
#  This is raw objective function type.
# 
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#  \item{expectedMean}{The expected mean (MxAlgebra object).}
#  \item{expectedCov}{The expected covariance matrix (MxAlgebra object).}
#  \item{observed}{The observed covariance matrix.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
#
#*/###########################################################################
setConstructorS3("RawObjective", function(expectedMean, expectedCov, observed) {

   if (missing(expectedMean)) expectedMean <- NA;
   if (missing(expectedCov)) expectedCov <- NA;
   if (missing(observed)) observed <- NA;

   extend(MxObjective(), "RawObjective",
	.expectedMean=expectedMean,
	.expectedCov=expectedCov,
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
setMethodS3("createMxJobClosureR", "RawObjective", function(this, job, ...) {

   model <- job$.model;
   model$makeFreeParametersList(); 

   function() {

      startValues <- model$getFreeParameters();
      expectedMean <- this$.expectedMean; # this is an MxAlgebra statement
      expectedCov <- this$.expectedCov;   # this is an MxAlgebra statement
      observed <- this$.observed;         # this is a matrix
      expectedMean$translateAlgebra();
      expectedCov$translateAlgebra();


      missdmnorm <- function(x, mu, sigma) {
         tsel <- !is.na(x);
         if (length(x[tsel]) == 0) return(NA);
         return(-2 * dmnorm(x[tsel], mu[tsel], sigma[tsel,tsel], log=TRUE));
      }

      objectiveFunction <- function(par) {
         model$updateMatrices(par);
         expectedMean$setDirtyBit(TRUE);
         expectedCov$setDirtyBit(TRUE);
         predictedMean <- expectedMean$evaluateTranslation();
         predictedCov <- expectedCov$evaluateTranslation();
         functionValue <- sum(apply(observed, 1, missdmnorm, mu=predictedMean, sigma=predictedCov), na.rm=TRUE);
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