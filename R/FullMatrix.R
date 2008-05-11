###########################################################################/**
# @RdocClass FullMatrix
#
# @title "The FullMatrix class"
#
# \description{
#
#  Creates a full matrix.
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
setConstructorS3("FullMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);

   valuesMatrix <- matrix(0, row, col);
   modifiable <- (row * col);
   freeParameters <- matrix(MxMatrix$FIXED(), row, col);
   if (free) freeParameters <- matrix(replicate(
      row * col, MxMatrix$createUniqueName()), row, col, byrow=TRUE);

   extend(MxMatrix(freeParameters, valuesMatrix, modifiable), "FullMatrix");

})

setMethodS3("setValuesWithList", "FullMatrix", function(this, valuesList,...) {
   this$.values <- t(matrix(valuesList, ncol(this$.values), nrow(this$.values)));
})

setMethodS3("setParametersWithList", "FullMatrix", 
   function(this, parametersList,...) {
     newMatrix <- matrix(parametersList, ncol(this$.parameters), nrow(this$.parameters));
     this$.parameters <- t(newMatrix);
})

setMethodS3("checkValidMatrix", "FullMatrix", 
   function(this, aMatrix,...) { TRUE })


setMethodS3("checkValidSpecification", "FullMatrix", 
   function(this, aMatrix,...) { TRUE })