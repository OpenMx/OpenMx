###########################################################################/**
# @RdocClass LowerMatrix
#
# @title "The LowerMatrix class"
#
# \description{
#
#  Creates a lower triangular square matrix.
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
setConstructorS3("LowerMatrix", function(row, col, free = FALSE) {

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
      ones[lower.tri(ones)] <- 0;
      # And set the remaining elements to 1, 2, ... , n
      ones[!lower.tri(ones)] <- 1:sum(ones);
      ones <- matrix(ones, row, col, byrow=TRUE);
      freeParameters <- ones;
   }

   extend(MxMatrix(freeParameters, valuesMatrix, modifiable), "LowerMatrix");

})


#########################################################################/**
# @RdocMethod checkValidMatrix
#
# @title "Inspect a Data Matrix for Validity"
#
# \seealso{
#     \link{MxMatrix.checkValidMatrix}
# }
#*/######################################################################### 
setMethodS3("checkValidMatrix", "LowerMatrix", function(this, aMatrix,...) {
   utriangle <- all(aMatrix[upper.tri(aMatrix, diag=FALSE)] == 0);
   return(utriangle);
})


#########################################################################/**
# @RdocMethod checkValidSpecification
#
# @title "Inspect a Parameter Matrix for Validity"
#
# \seealso{
#     \link{MxMatrix.checkValidSpecification}
# }
#*/######################################################################### 
setMethodS3("checkValidSpecification", "LowerMatrix", function(this, aMatrix,...) {
   utriangle <- all(aMatrix[upper.tri(aMatrix, diag=FALSE)] == 0);
   return(utriangle);
})


#########################################################################/**
# @RdocMethod setValuesWithList
#
# @title "Set Matrix Values With A List"
#
# \seealso{
#     \link{MxMatrix.setValuesWithList}
# }
#*/######################################################################### 
setMethodS3("setValuesWithList", "LowerMatrix", function(this, valuesList,...) {
   # Set the upper triangle to valuesList
   this$values[upper.tri(this$values, diag=TRUE)] <- valuesList;
   # And now perform a transpose operation
   this$values <- t(this$values);
   this$values[upper.tri(this$values, diag=FALSE)] <- 0;
})


#########################################################################/**
# @RdocMethod setParametersWithList
#
# @title "Set Matrix Parameters With A List"
#
# \seealso{
#     \link{MxMatrix.setParametersWithList}
# }
#*/######################################################################### 
setMethodS3("setParametersWithList", "LowerMatrix", function(this, parametersList,...) {
   # Set the upper triangle to parametersList
   this$parameters[upper.tri(this$parameters, diag=TRUE)] <- parametersList;
   # And now perform a transpose operation
   this$parameters <- t(this$parameters);
   this$parameters[upper.tri(this$parameters, diag=FALSE)] <- 0;
})
