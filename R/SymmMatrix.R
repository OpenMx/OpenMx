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
      # Set the lower triangular matrix to zero
      ones[lower.tri(ones)] <- 0;
      # And set the remaining elements to 1, 2, ... , n
      uniqueNames <- replicate(sum(ones), MxMatrix$createUniqueName());
      target <- matrix("", row, col); 
      # Set the upper triangle to uniqueNames
      target[upper.tri(target, diag=TRUE)] <- uniqueNames;
      # And now perform a transpose operation
      target <- t(target);
   
      triangle <- matrix("", nrow(target), ncol(target));
      triangle[lower.tri(triangle)] <- target[lower.tri(target)];

      pasteParameters <- mapply(function(x,y) {paste(x,y,sep="")},
         target, t(triangle));

      freeParameters <- matrix(pasteParameters, 
         row, col);
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
   this$.values[upper.tri(this$.values, diag=TRUE)] <- valuesList;
   this$.values <- this$.values + t(this$.values) - diag(diag(this$.values));
})


setMethodS3("setParametersWithList", "SymmMatrix", function(this, parametersList,...) {
   target <- matrix("", nrow(this$.parameters), ncol(this$.parameters)); 
   # Set the upper triangle to parametersList
   target[upper.tri(target, diag=TRUE)] <- parametersList;
   # And now perform a transpose operation
   target <- t(target);
   
   triangle <- matrix("", nrow(target), ncol(target));
   triangle[lower.tri(triangle)] <- target[lower.tri(target)];

   pasteParameters <- mapply(function(x,y) {paste(x,y,sep="")},
      target, t(triangle));

   this$.parameters <- matrix(pasteParameters, 
      nrow(this$.parameters), ncol(this$.parameters));

})
