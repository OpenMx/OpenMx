###########################################################################/**
# @RdocClass StandMatrix
#
# @title "The StandMatrix class"
#
# \description{
#
#  Creates a standardized matrix. This is a symmetric square matrix
#  with ones along the diagonal.
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
setConstructorS3("StandMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkSquare(row,col);

   valuesMatrix <- diag(row);
   triangle <- lower.tri(valuesMatrix, diag = FALSE);
   modifiable <- length(triangle[triangle]);
   freeParameters <- matrix(MxMatrix$FIXED(), row, col);
   if (free) {
      ones <- matrix(1, row, col);
      # Set the lower triangular matrix to zero
      ones[lower.tri(ones, diag=TRUE)] <- 0;
      # And set the remaining elements to 1, 2, ... , n
      uniqueNames <- replicate(sum(ones), MxMatrix$createUniqueName());
      freeParameters <- matrix("", nrow(freeParameters), ncol(freeParameters)); 
      # Set the upper triangle to parametersList
      freeParameters[upper.tri(freeParameters)] <- uniqueNames;
      # And now perform a transpose operation
      freeParameters <- t(freeParameters);
      pasteParameters <- mapply(function(x,y) {paste(x,y,sep="")},
         freeParameters, t(freeParameters))
      freeParameters <- matrix(pasteParameters, 
         nrow(freeParameters), ncol(freeParameters))
      diag(freeParameters) <- MxMatrix$FIXED();
   }

   extend(MxMatrix(freeParameters, valuesMatrix, modifiable), "StandMatrix");

})


setMethodS3("checkValidMatrix", "StandMatrix", function(this, aMatrix,...) {
   symmetry <- all(aMatrix == t(aMatrix));
   ones <- all(diag(aMatrix) == 1);
   return(symmetry && ones);
})


setMethodS3("checkValidSpecification", "StandMatrix", function(this, aMatrix,...) {
   symmetry <- all(aMatrix == t(aMatrix));
   ones <- all(diag(aMatrix) == MxMatrix$FIXED());
   return(symmetry && ones);
})


setMethodS3("setValuesWithList", "StandMatrix", function(this, valuesList,...) {

   # Set the lower triangular matrix to zero
   this$.values[lower.tri(this$.values, diag=TRUE)] <- 0;

   # And set the remaining elements to valuesList
   this$.values[!lower.tri(this$.values, diag=TRUE)] <- valuesList;
   this$.values <- matrix(this$.values, nrow(this$.values), ncol(this$.values), byrow=TRUE);
   this$.values <- this$.values + t(this$.values) + diag(nrow = nrow(this$.values));
})


setMethodS3("setParametersWithList", "StandMatrix", function(this, parametersList,...) {
   target <- matrix("", nrow(this$.parameters), ncol(this$.parameters)); 
   # Set the upper triangle to parametersList
   target[upper.tri(target)] <- parametersList;

   pasteParameters <- mapply(function(x,y) {paste(x,y,sep="")},
      target, t(target))

   this$.parameters <- matrix(pasteParameters, 
      nrow(this$.parameters), ncol(this$.parameters))
   diag(this$.parameters) <- MxMatrix$FIXED();

})
