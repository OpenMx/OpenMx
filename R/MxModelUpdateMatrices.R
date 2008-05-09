
   # DO NOT export this method
   # this is a helper function to updateMatrices()
   updateMatrix <- function(reference, row, col, newvalue) {
      reference$values[row,col] <- newvalue;
   };

   # DO NOT export this method
   # this is a helper function to updateMatrices()
   updateMatricesHelper <- function(listTuples, parameter) {
      lapply(listTuples, function(aTuple) {
         updateMatrix(aTuple[[1]], aTuple[[2]], aTuple[[3]], parameter)
      });
   };

   # DO NOT export this method
   # length(parameters) must be equal to length(freeVariables)
   setMethodS3("updateMatrices", "MxModel", function(this, parameters, ...) {
      returnValue <- mapply(updateMatricesHelper, this$.freeVariablesList, parameters);
      invisible(returnValue);
   });
