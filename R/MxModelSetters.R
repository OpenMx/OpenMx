
setMethodS3("transformMatrix", "MxMatrix", function(this, offset,...) {
  newFree <- this$parameters
  newFree[newFree > 0] <- newFree[newFree > 0] + offset
  this$parameters <- newFree
})


removeNullElements <- function(s) {
  return(s[!sapply(s, is.null)])   
}

setMethodS3("updateFreeVariablesList", "MxModel", function(this,...) {
  updatedList <- list();
  for (aField in this$getFields()) {
    if (inherits(this[[aField]],"MxMatrix")) {
      parameters <- this[[aField]]$parameters;
      for(row in 1:dim(parameters)[1]) {
        for(col in 1:dim(parameters)[2]) {
          variableNumber <- parameters[row,col];
          if (variableNumber > 0) {
            if (variableNumber > length(updatedList) || 
              is.null(updatedList[[variableNumber]])) {
              updatedList[[variableNumber]] <- list(list(matrix=this[[aField]], row=row, col=col));
            } else {
              copies <- length(updatedList[[variableNumber]]);
              updatedList[[variableNumber]][[copies + 1]] <- list(matrix=this[[aField]], row=row, col=col);
            }
          }
        }
      }
    }
  }
  updatedList <- removeNullElements(updatedList);
  this$.freeVariablesList <- updatedList;
})
 

setMethodS3("setValues", "MxMatrix", function(this, values,...) {

   if (is.vector(values)) {
      modifiable <- this$.modifiable;
      if (modifiable != length(values)) {
         error <- paste("Your values list has", length(values),
            	"elements but matrix has", modifiable,
            	"modifiable elements.");
         throw(error);
      }
      this$setValuesWithList(values);
   } else if (is.matrix(values)) {
      if (!all(dim(this$values) == dim(values))) {
          error <- paste("Second argument has dimensions",
            	paste(dim(values), collapse = " "), "but matrix",
            	"has dimensions",
            	paste(dim(this$values), collapse = " "), ".");
          throw(error);
      }
      valid <- this$checkValidMatrix(values);
      if (!valid) {
          error <- paste("Second argument is not a valid",
            	data.class(this),".");
          throw(error);
      }
      this$values <- values;
   } else {
      throw("Second argument is neither a vector nor a matrix.");
   }
})

setMethodS3("setParameters", "MxMatrix", function(this, parameters,...) {

   if (is.vector(parameters)) {
      modifiable <- this$.modifiable;
      if (modifiable != length(parameters)) {
         error <- paste("Your values list has", length(parameters),
            	"elements but matrix has", modifiable,
            	"modifiable elements.");
         throw(error);
      }
      this$setParametersWithList(parameters);
   } else if (is.matrix(parameters)) {
      if (!all(dim(this$parameters) == dim(parameters))) {
          error <- paste("Second argument has dimensions",
            	paste(dim(parameters), collapse = " "), "but matrix",
            	"has dimensions",
            	paste(dim(matrix$parameters), collapse = " "), ".")
          throw(error);
      }
      valid <- this$checkValidMatrix(parameters);
      if (!valid) {
          error <- paste("Second argument is not a valid",
            	data.class(this),".");
          throw(error);
      }
      this$parameters <- parameters
   } else {
      throw("Second argument is neither a vector nor a matrix.");
   }

})
