###########################################################################/**
# @RdocClass IdenMatrix
#
# @title "The IdenMatrix class"
#
# \description{
#
#  Creates a square identity matrix. This matrix has no
#  free parameters.
#
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{row}{The number of rows this matrix contains.}
#   \item{col}{The number of columns this matrix contains.}
# }
#
#
#
#*/###########################################################################
setConstructorS3("IdenMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkFree(free);
   checkSquare(row,col);

   freeParameters <- matrix(0, row, col);
   valuesMatrix <- diag(row);
   modifiable <- 0;
   extend(MxMatrix(freeParameters, valuesMatrix, modifiable), "IdenMatrix");

})
