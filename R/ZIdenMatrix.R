###########################################################################/**
# @RdocClass ZIdenMatrix
#
# @title "The ZIdenMatrix class"
#
# \description{
#
#  Creates a partitioned zero|identity matrix. This matrix has no
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
setConstructorS3("ZIdenMatrix", function(row, col) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);

   freeParameters <- matrix(MxMatrix$FIXED(), row, col);
   if (row == col) {
      valuesMatrix <- diag(row);
   } else if (row < col) {
      square <- diag(row);
      zeros  <- matrix(0, row, col - row);
      valuesMatrix <- cbind(zeros, square);
   } else {
      square <- diag(col);
      zeros  <- matrix(0, col, row - col);
      valuesMatrix <- rbind(zeros, square);
   }
   modifiable <- 0;
   extend(MxMatrix(freeParameters, valuesMatrix, modifiable), "ZIdenMatrix");

})
