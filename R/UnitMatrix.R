###########################################################################/**
# @RdocClass UnitMatrix
#
# @title "The UnitMatrix class"
#
# \description{
#
#  Creates a matrix that contains only unity values. This matrix has no
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
setConstructorS3("UnitMatrix", function(row, col) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);

   freeParameters <- matrix(0, row, col);
   valuesMatrix <- matrix(1, row, col);
   modifiable <- 0;
   extend(MxMatrix(freeParameters, valuesMatrix, modifiable), "UnitMatrix");

})
