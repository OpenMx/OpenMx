setConstructorS3("MxAlgebra", function(formula) {

  if (missing(formula)) formula <- NA;
  formula <- match.call()$formula;

  extend(Object(), "MxAlgebra",
	.formula=formula,
	.translation=NULL,
	.dirty=TRUE,
	.value=0
  );

})

setMethodS3("print", "MxAlgebra", function(x, ...) {
   cat("Formula: ");
   print(x$.formula);
   cat("Translation: ");
   print(x$.translation);
   cat("Dirty:", x$.dirty, "\n");
   if (!x$.dirty) {
      cat("Value:" );
      print(x$.value);
   }
   invisible(x)
})

translateAlgebraHelper <- function(formula) {
   if(is.call(formula) && formula[[1]] != '$') {
      return(as.call(lapply(formula, translateAlgebraHelper)));
   }
   evalFormula <- eval(formula);
   if (is.object(evalFormula)) {
      if (inherits(evalFormula,"MxAlgebra")) {
         if (is.null(evalFormula$.translation)) {
            evalFormula$translateAlgebra();
         }
         return(substitute(x$evalTranslation(), list(x = formula)));
      return();
      } else if (inherits(evalFormula,"MxMatrix")) {
         return(substitute(x$values, list(x = formula)));
      }
   }
   return(formula);
};

setMethodS3("translateAlgebra", "MxAlgebra", function(this, ...) {
   this$.translation <- translateAlgebraHelper(this$.formula);
});


setMethodS3("evalTranslation", "MxAlgebra", function(this, ...) {
   if(!this$.dirty) {
      return(this$.value);
   }
   value <- eval(this$.translation);
   this$.value <- value;
   this$.dirty <- FALSE;
   return(value);
});

setDirtyBitHelper <- function(formula, value) {
   if(is.call(formula) && formula[[1]] != '$') {
      lapply(formula, setDirtyBitHelper);
   }
   evalFormula <- eval(formula);
   if (is.object(evalFormula)) {
      if (inherits(evalFormula,"MxAlgebra")) {
         evalFormula$.dirty <- value;
      }
   }
};

setMethodS3("setDirtyBit", "MxAlgebra", function(this, value, ...) {
   this$.dirty <- value;
   setDirtyBitHelper(this$.translation, value);
});
