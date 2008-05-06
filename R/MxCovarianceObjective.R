setConstructorS3("CovarianceObjective", function(expected, observed) {

   if (missing(expected)) row <- NA;
   if (missing(observed)) col <- NA;


   extend(MxObjective(), "CovarianceObjective",
	.expected=expected,
	.observed=observed);

})


setMethodS3("createMxJobClosureR", "MxCovarianceObjective", function(objective, job, ...) {
      
   model <- job$.model; 

   function() {

      observedCov <- objective$.observed;
      startValues <- model$getFreeVariables();

      objectiveFunction <- function(par) {
         model$updateMatrices(par);
         predictedCov <- model$evaluateAlgebra(objective$.expected);
         functionValue <- sum(diag(observedCov %*% solve(predictedCov))) + log(det(predictedCov));
         functionValue;
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