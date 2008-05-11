###########################################################################/**
# @RdocClass MxModel
#
# @title "The MxModel class"
#
# \description{
#
#  Foo.
# 
#  @classhierarchy
# }
# 
# @synopsis
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
#
#*/###########################################################################
setConstructorS3("MxModel", function() { 
  
  freeParametersList <- list();
  
  extend(Object(), "MxModel", 
      .freeParametersList = freeParametersList);

})




###########################################################################/**
# @RdocFunction updateMatricesHelper
#
# @title "UpdateMatricesHelper"
# \description{
#   This is a helper function to MxModel.updateMatrices.
# }
#
#*/######################################################################### 
updateMatricesHelper <- function(listTuples, parameter) {
   lapply(listTuples, function(aTuple) {
      updateMatrix(aTuple[[1]], aTuple[[2]], aTuple[[3]], parameter)
   });
};


#########################################################################/**
# @RdocMethod updateMatrices
#
# @title "Update Values In Objective Function"
# 
# \description{
#    Given a vector of new values, update the free parameters of the
#    model with these values.  The size of the parameters vector must
#    equal the size of the free parameters list. 
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxModel object.}
#  \item{parameters}{A vector of new values.}
#  \item{...}{Unused.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("updateMatrices", "MxModel", function(this, parameters, ...) {
   returnValue <- mapply(updateMatricesHelper, 
      this$.freeParametersList, parameters);
   invisible(returnValue);
});


# DO NOT export this method
# this is a helper function to getFreeParameters()
getFreeParametersMatrix <- function(threeTuple) {
   reference <- threeTuple[[1]];
   row <- threeTuple[[2]];
   col <- threeTuple[[3]];
   return(reference$.values[row,col]);
}

# DO NOT export this method
# this is a helper function to getFreeParameters()
#
# Only bother getting the first value,
# since multiple matrix locations for the same free parameter
# should have identical values.
getFreeParametersHelper <- function(listTuples) {
   return(getFreeParametersMatrix(listTuples[[1]]));
}

#
# EXPORT this method
#
setMethodS3("getFreeParameters", "MxModel", function(this,...) {
   return(unlist(lapply(this$.freeParametersList, getFreeParametersHelper)))
});



setMethodS3("makeFreeParametersList", "MxModel", function(this,...) {
  freeParameters <- list();
  for (aField in this$getFields()) {
    if (inherits(this[[aField]],"MxMatrix")) {
      pMatrix <- this[[aField]]$.parameters;
      for(row in 1:dim(pMatrix)[1]) {
        for(col in 1:dim(pMatrix)[2]) {
          parameter <- pMatrix[row,col];
          if (parameter != MxMatrix$FIXED()) {
            if (is.null(freeParameters[[parameter]])) {
              freeParameters[[parameter]] <- list(list(matrix=this[[aField]], row=row, col=col));
            } else {
              size <- length(freeParameters[[parameter]]);
              freeParameters[[parameter]][[size + 1]] <- list(matrix=this[[aField]], row=row, col=col);
            }
          }
        }
      }
    }
  }
  this$.freeParametersList <- freeParameters;
  return(freeParameters);
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
      if (!all(dim(this$.values) == dim(values))) {
          error <- paste("Second argument has dimensions",
            	paste(dim(values), collapse = " "), "but matrix",
            	"has dimensions",
            	paste(dim(this$.values), collapse = " "), ".");
          throw(error);
      }
      valid <- this$checkValidMatrix(values);
      if (!valid) {
          error <- paste("Second argument is not a valid",
            	data.class(this),".");
          throw(error);
      }
      this$.values <- values;
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
      this$setParametersWithList(as.character(parameters));
   } else if (is.matrix(parameters)) {
      if (!all(dim(this$.parameters) == dim(parameters))) {
          error <- paste("Second argument has dimensions",
            	paste(dim(parameters), collapse = " "), "but matrix",
            	"has dimensions",
            	paste(dim(this$.parameters), collapse = " "), ".")
          throw(error);
      }
      valid <- this$checkValidMatrix(parameters);
      if (!valid) {
          error <- paste("Second argument is not a valid",
            	data.class(this),".");
          throw(error);
      }
      this$.parameters <- matrix(as.character(parameters), nrow(parameters), ncol(parameters))
   } else {
      throw("Second argument is neither a vector nor a matrix.");
   }

})
