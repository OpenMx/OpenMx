###########################################################################/**
# @RdocClass MxMatrix
#
# @title "The MxMatrix class"
#
# \description{
#
#  This is the abstract superclass of all MxMatrix types.
# 
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{parameters}{A @matrix of parameter strings.}
#   \item{values}{A @matrix of initial values.}
#   \item{modifiable}{The number of modifiable elements in this matrix type.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
#
#*/###########################################################################
setConstructorS3("MxMatrix", function(parameters, values, modifiable) {

  if (missing(parameters)) parameters <- NA;
  if (missing(values))  values <- NA;
  if (missing(modifiable)) modifiable <- NA;

  extend(Object(), "MxMatrix",
    parameters=parameters,
    values=values,
    .modifiable=modifiable
  );

}, abstract=TRUE)

#########################################################################/**
# @RdocMethod print
#
# @title "Print MxMatrix"
# 
# \description{
#   \code{print} prints its matrix and returns
#      it invisibly (via \code{invisible}(x)).
# }
#
# @synopsis
#
# \arguments{
#  \item{x}{The MxMatrix object.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("print", "MxMatrix", function(x, ...) {
   cat(paste("MxMatrix:", data.class(x)), sep="\n")
   cat("Parameters: ", sep="\n")
   print(x$parameters)
   cat("Values: ", sep="\n")
   print(x$values)
   cat("Modifiable:", x$.modifiable)
   cat("\n")
   invisible(x)
})


checkMatrix <- function(row,col) {
   if (row < 0) throw("Row value is a negative number");
   if (col < 0) throw("Col value is a negative number");
   if (trunc(row) != row) throw("Row value is not an integer");
   if (trunc(col) != col) throw("Col value is not an integer");
}


checkSquare <- function(row, col) {
   if (row != col) throw("Row and column dimensions do not match");
}


#TODO IMPLEMENT
#SDiagMatrix is a subdiagonal matrix (zeros on and above the diagonal)
#
