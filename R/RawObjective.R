setConstructorS3("RawObjective", function(expectedMean, expectedCov, observed) {

   if (missing(expectedMean)) expectedMean <- NA;
   if (missing(expectedCov)) expectedCov <- NA;
   if (missing(observed)) observed <- NA;

   extend(MxObjective(), "RawObjective",
	.expectedMean=expectedMean,
	.expectedCov=expectedCov,
	.observed=observed);

})

setMethodS3("createMxJobClosureR", "RawObjective", function(objective, job, ...) {

   model <- job$.model;
   model$makeFreeParametersList(); 

   function() {

      startValues <- model$getFreeParameters();
      expectedMean <- objective$.expectedMean; # this is an MxAlgebra statement
      expectedCov <- objective$.expectedCov;   # this is an MxAlgebra statement
      observed <- objective$.observed;         # this is a matrix
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
         predictedMean <- expectedMean$evalTranslation();
         predictedCov <- expectedCov$evalTranslation();
         functionValue <- sum(apply(observed, 1, missdmnorm, mu=predictedMean,
		sigma=predictedCov), na.rm=TRUE);
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