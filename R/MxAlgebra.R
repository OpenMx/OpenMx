setConstructorS3("MxAlgebra", function(formula) {

  if (missing(formula)) formula <- NA;
  formula <- match.call()$formula;

  extend(Object(), "MxAlgebra",
	.formula=formula,
	.translation=NULL,
	.dirty=TRUE
  );

})

setMethodS3("print", "MxAlgebra", function(x, ...) {
   cat("Formula: ")
   print(x$.formula)
   cat("Translation: ")
   print(x$.translation)
   cat("Dirty:", x$.dirty, "\n")
   invisible(x)
})
