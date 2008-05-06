#
# MxJob is the abstract supertype of all job types
#
setConstructorS3("MxJob", function(model, objective) {

   if (missing(model)) model <- NA;
   if (missing(objective)) objective <- NA;


   extend(Object(), "MxJob",
	.model=model,
	.objective=objective);

})


createMxClosure <- function(job, use_R = FALSE) {
   if (use_R) {
      createMxJobClosureR(job$.objective, job);
   } else {
      createMxJobClosureC(job$.objective, job);
   }
}

