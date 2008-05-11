###########################################################################/**
# @RdocClass DiagMatrix
#
# @title "The DiagMatrix class"
#
# \description{
#
#  Creates a square diagonal matrix. Parameters and values along the
#  main diagonal can be modified.
#
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{row}{The number of rows this matrix contains.}
#   \item{col}{The number of columns this matrix contains.}
#   \item{free}{TRUE if this matrix has free parameters.}
# }
#
#
#
#*/###########################################################################
setConstructorS3("DiagMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkSquare(row,col);

   freeParameters <- matrix(0, row, col);
   if (free) freeParameters <- diag(row) * 1:row;
   valuesMatrix <- matrix(0, row, col);
   modifiable <- row;

   extend(MxMatrix(freeParameters, valuesMatrix, modifiable), "DiagMatrix");

})

setMethodS3("checkValidMatrix", "DiagMatrix", function(this, aMatrix,...) {
   ltriangle <- all(aMatrix[lower.tri(aMatrix, diag=FALSE)] == 0);
   utriangle <- all(aMatrix[upper.tri(aMatrix, diag=FALSE)] == 0);
   return(ltriangle && utriangle);
})

setMethodS3("checkValidSpecification", "DiagMatrix", function(this, aMatrix,...) {
   ltriangle <- all(aMatrix[lower.tri(aMatrix, diag=FALSE)] == 0);
   utriangle <- all(aMatrix[upper.tri(aMatrix, diag=FALSE)] == 0);
   return(ltriangle && utriangle);
})

setMethodS3("setValuesWithList", "DiagMatrix", function(this, valuesList,...) {
   this$values <- diag(valuesList);
})

setMethodS3("setParametersWithList", "DiagMatrix", function(this, parametersList,...) {
   this$parameters <- diag(parametersList);
})
