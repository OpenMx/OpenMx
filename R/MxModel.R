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

#########################################################################/**
# @RdocFunction getFreeParameterValue
#
# @title "Helper Function to getFreeParametersHelper"
# 
# \description{    
#    Returns the current value of the free parameter for the 
#    (MxMatrix, row, col) specification passed an argument.
# }
#
# @synopsis
#
# \arguments{
#  \item{threeTuple}{One 3-tuples that consists of(MxMatrix, row, column).}
# }
#
#*/######################################################################### 
getFreeParameterValue <- function(threeTuple) {
   reference <- threeTuple[[1]];
   row <- threeTuple[[2]];
   col <- threeTuple[[3]];
   return(reference$.values[row,col]);
}

#########################################################################/**
# @RdocFunction getFreeParametersHelper
#
# @title "Helper Function to getFreeParameters"
# 
# \description{
#    
#    Only bother getting the first value,
#    since multiple matrix locations for the same free parameter
#    should have identical values.
# }
#
# @synopsis
#
# \arguments{
#  \item{listTuples}{A list of 3-tuples (MxMatrix, row, column).}
# }
#
#*/######################################################################### 
getFreeParametersHelper <- function(listTuples) {
   return(getFreeParameterValue(listTuples[[1]]));
}

#########################################################################/**
# @RdocMethod getFreeParameters
#
# @title "Return a Vector of Free Parameter Values"
# 
# \description{
#    Returns the current values of the model's free parameters.
#    This function uses the \code{.freeParametersList} field, which
#    requires that the \link{MxModel.makeFreeParametersList} method
#    has been called before this method is called.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxModel object.}
#  \item{...}{Unused.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getFreeParameters", "MxModel", function(this,...) {
   return(unlist(lapply(this$.freeParametersList, getFreeParametersHelper)))
});


#########################################################################/**
# @RdocMethod makeFreeParametersList
#
# @title "Make the Free Parameter List"
# 
# \description{
#    Constructs the free parameter list, stores the list in the 
#    \code{.freeParametersList} field, and returns a copy of this list.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxModel object.}
#  \item{...}{Unused.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
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
 

