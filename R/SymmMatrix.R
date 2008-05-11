###########################################################################/**
# @RdocClass SymmMatrix
#
# @title "The SymmMatrix class"
#
# \description{
#
#  Creates a symmetric square matrix.
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
setConstructorS3("SymmMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkSquare(row,col);

   valuesMatrix <- matrix(0, row, col);
   triangle <- lower.tri(valuesMatrix, diag = TRUE);
   modifiable <- length(triangle[triangle]);
   freeParameters <- matrix(0, row, col);
   if (free) {
      ones <- matrix(1, row, col);
      # Set the upper triangular matrix to zero
      ones[lower.tri(ones, diag=TRUE)] <- 0;
      # And set the remaining elements to 1, 2, ... , n
      ones[!lower.tri(ones, diag=TRUE)] <- 1:sum(ones);
      # And now fill in the diagonals
      biggest <- max(ones);
      diagonals <- diag(row) * (biggest + 1) : (biggest + row);
      freeParameters <- ones + t(ones) + diagonals;
   }

   extend(MxMatrix(freeParameters, valuesMatrix, modifiable), "SymmMatrix");

})

setMethodS3("checkValidMatrix", "SymmMatrix", function(this, aMatrix,...) {
   symmetry <- all(aMatrix == t(aMatrix));
   return(symmetry);
})


setMethodS3("checkValidSpecification", "SymmMatrix", function(this, aMatrix,...) {
   symmetry <- all(aMatrix == t(aMatrix));
   return(symmetry);
})


setMethodS3("setValuesWithList", "SymmMatrix", function(this, valuesList,...) {
   # Set the lower triangle to valuesList
   this$values[upper.tri(this$values, diag=TRUE)] <- valuesList;
   this$values <- this$values + t(this$values) - diag(diag(this$values));
})


setMethodS3("setParametersWithList", "SymmMatrix", function(this, parametersList,...) {
   # Set the lower triangle to parametersList
   this$parameters[upper.tri(this$parameters, diag=TRUE)] <- parametersList;
   this$parameters <- this$parameters + t(this$parameters) - diag(diag(this$parameters));
})
