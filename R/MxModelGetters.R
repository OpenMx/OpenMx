
   # DO NOT export this method
   # this is a helper function to getFreeVariables()
   getFreeVariablesMatrix <- function(threeTuple) {
      reference <- threeTuple[[1]];
      row <- threeTuple[[2]];
      col <- threeTuple[[3]];
      return(reference$values[row,col]);
   }

   # DO NOT export this method
   # this is a helper function to getFreeVariables()
   #
   # Only bother getting the first value,
   # since multiple matrix locations for the same free parameter
   # should have identical values.
   getFreeVariablesHelper <- function(listTuples) {
      return(getFreeVariablesMatrix(listTuples[[1]]));
   }

   #
   # EXPORT this method
   #
   setMethodS3("getFreeVariables", "MxModel", function(this,...) {
      return(unlist(lapply(this$.freeVariablesList, getFreeVariablesHelper)))
   });