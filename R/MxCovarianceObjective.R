setConstructorS3("CovarianceObjective", function(expected, observed) {

   if (missing(expected)) expected <- NA;
   if (missing(observed)) observed <- NA;


   extend(MxObjective(), "CovarianceObjective",
	.expected=expected,
	.observed=observed);

})


setMethodS3("createMxJobClosureR", "CovarianceObjective", function(objective, job, ...) {
      
   model <- job$.model;
   model$updateFreeVariablesList(); 

   function() {

      observedCov <- objective$.observed; # this is a matrix
      expectedCov <- objective$.expected; # this is an MxAlgebra statement
      startValues <- model$getFreeVariables();
      expectedCov$translateAlgebra();

      objectiveFunction <- function(par) {
         model$updateMatrices(par);
         expectedCov$setDirtyBit(TRUE);
         predictedCov <- expectedCov$evalTranslation();
         functionValue <- sum(diag(observedCov %*% solve(predictedCov))) + log(det(predictedCov));
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
