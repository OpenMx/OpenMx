setConstructorS3("MxAlgebra", function(formula) {

  if (missing(formula)) formula <- NA;
  formula <- match.call()$formula;

  extend(Object(), "MxAlgebra",
	.formula=formula,
	.translation=NULL,
	.dirty=TRUE
  );

})

setMethodS3("print", "MxAlgebra", function(this, ...) {
   cat("Formula: ")
   print(this$.formula)
   cat("Translation: ")
   print(this$.translation)
   cat("Dirty:", this$.dirty, "\n")
   invisible(this)
})
