###########################################################################/**
# @RdocClass MxAlgebra
#
# @title "The MxAlgebra class"
#
# \description{
#
#  This class stores a matrix algebra expression to be evaluated
#  at some later point.
# 
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{formula}{A matrix algebra expression that
#                  will be evaluated at some later point.}
# }
#
# \section{Fields and Methods}{
#  \bold{Fields:}
#  \itemize{
#     \item \code{.formula} - The untranslated matrix algebra expression.
#     This expression cannot be evaluated until it has been translated.
#     \item \code{.translation} - The translated matrix algebra expression.
#     This field is populated by a call to 
#     \code{\link[=MxAlgebra.translateAlgebra]{translateAlgebra}}
#     and it is evaluated by a call to 
#     \code{\link[=MxAlgebra.evaluateTranslation]{evaluateTranslation}}.
#     \item \code{.dirty} - If \code{TRUE},
#     then \code{\link[=MxAlgebra.evaluateTranslation]{evaluateTranslation}}
#     will used the cached value in \code{.value}. This field must be set
#     explicitly to \code{TRUE} or \code{FALSE} by a call to
#     \code{\link[=MxAlgebra.setDirtyBit]{setDirtyBit}}.  
#     \item \code{.value} - The cached value of the MxAlgebra expression.
#  }
#  @allmethods
# }
# 
#
#*/###########################################################################
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

#########################################################################/**
# @RdocMethod print
#
# @title "Print MxAlgebra"
# 
# \description{
#   \code{print} prints the matrix algebra object and returns
#      it invisibly (via \code{invisible}(x)).
# }
#
# @synopsis
#
# \arguments{
#  \item{x}{The MxAlgebra object.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
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


#########################################################################/**
# @RdocFunction translateAlgebraHelper
#
# @title "Helper Function to translateAlgebra"
# 
# \description{
#    A helper function that recursively translates all MxAlgebra 
#    objects contained within the formula argument.
# }
#
# @synopsis
#
# \arguments{
#  \item{formula}{The formula to search for MxAlgebras.}
# }
#
#*/######################################################################### 
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
         return(substitute(x$evaluateTranslation(), list(x = formula)));
      return();
      } else if (inherits(evalFormula,"MxMatrix")) {
         return(substitute(x$.values, list(x = formula)));
      }
   }
   return(formula);
};

#########################################################################/**
# @RdocMethod translateAlgebra
#
# @title "Translate MxAlgebra for Evaluation"
# 
# \description{
#    Translates an MxAlgebra object so that it can be evaluated
#    by a call to \link{MxAlgebra.evaluateTranslation}.
#    Also recursively translates all MxAlgebra objects contained
#    within the \code{.formula} field of this MxAlgebra object.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxAlgebra object.}
#  \item{...}{Unused.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("translateAlgebra", "MxAlgebra", function(this, ...) {
   this$.translation <- translateAlgebraHelper(this$.formula);
});


#########################################################################/**
# @RdocMethod evaluateTranslation
#
# @title "Evaluate a Translated MxAlgebra Expression"
# 
# \description{
#    Evaluates an MxAlgebra expression that has been translated
#    by a call to \link{MxAlgebra.translateAlgebra}.
#    If the dirty bit on this object is set to TRUE, then we have
#    to actually evaluate the expression.  If the dirty bit is set
#    to FALSE, then we just return the currently cached value.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxAlgebra object.}
#  \item{...}{Unused.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("evaluateTranslation", "MxAlgebra", function(this, ...) {
   if(!this$.dirty) {
      return(this$.value);
   }
   value <- eval(this$.translation);
   this$.value <- value;
   this$.dirty <- FALSE;
   return(value);
});

#########################################################################/**
# @RdocFunction setDirtyBitHelper
#
# @title "Helper Function to setDirtyBit"
# 
# \description{
#    A helper function that recursively sets the dirty bit 
#    to a new value for all MxAlgebra objects contained
#    within the formula argument.
# }
#
# @synopsis
#
# \arguments{
#  \item{formula}{The formula to search for MxAlgebras.}
#  \item{value}{The new boolean value.}
# }
#
#*/######################################################################### 
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

#########################################################################/**
# @RdocMethod setDirtyBit
#
# @title "Set Dirty Bit in MxAlgebra"
# 
# \description{
#    Sets the dirty bit in a MxAlgebra object to a new value.
#    Also recursively updates the dirty bit in all MxAlgebra objects contained
#    within the \code{.formula} field of this MxAlgebra object.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxAlgebra object.}
#  \item{value}{The new boolean value.  Must be either TRUE or FALSE.}
#  \item{...}{Unused.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setDirtyBit", "MxAlgebra", function(this, value, ...) {
   this$.dirty <- value;
   setDirtyBitHelper(this$.translation, value);
});
