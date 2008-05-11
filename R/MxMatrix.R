###########################################################################/**
# @RdocClass MxMatrix
#
# @title "The MxMatrix class"
#
# \description{
#
#  This is the abstract superclass of all MxMatrix types.
# 
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{parameters}{A @matrix of parameter strings.}
#   \item{values}{A @matrix of initial values.}
#   \item{modifiable}{The number of modifiable elements in this matrix type.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
#
#*/###########################################################################
setConstructorS3("MxMatrix", function(parameters, values, modifiable) {

  if (missing(parameters)) parameters <- NA;
  if (missing(values))  values <- NA;
  if (missing(modifiable)) modifiable <- NA;

  extend(Object(), "MxMatrix",
    .parameters=parameters,
    .values=values,
    .modifiable=modifiable,
    .uniqueCount=1    
  );

}, abstract = TRUE);


setMethodS3("FIXED", "MxMatrix", function(this, ...) {
   return("0");
}, static = TRUE);

setMethodS3("UNIQUE", "MxMatrix", function(this, ...) {
   return("auto-");
}, static = TRUE);


setMethodS3("createUniqueName", "MxMatrix", function(this, ...) {
   retval <- paste(MxMatrix$UNIQUE(), MxMatrix$.uniqueCount, sep="");
   MxMatrix$.uniqueCount <- MxMatrix$.uniqueCount + 1; 
   return(retval);
}); 


#########################################################################/**
# @RdocMethod print
#
# @title "Print MxMatrix"
# 
# \description{
#   \code{print} prints its matrix and returns
#      it invisibly (via \code{invisible}(x)).
# }
#
# @synopsis
#
# \arguments{
#  \item{x}{The MxMatrix object.}
# }
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("print", "MxMatrix", function(x, ...) {
   cat(paste("MxMatrix:", data.class(x)), sep="\n")
   cat("Parameters: ", sep="\n")
   print(x$.parameters)
   cat("Values: ", sep="\n")
   print(x$.values)
   cat("Modifiable:", x$.modifiable)
   cat("\n")
   invisible(x)
})


checkMatrix <- function(row,col) {
   if (row < 0) throw("Row value is a negative number");
   if (col < 0) throw("Col value is a negative number");
   if (trunc(row) != row) throw("Row value is not an integer");
   if (trunc(col) != col) throw("Col value is not an integer");
}


checkSquare <- function(row, col) {
   if (row != col) throw("Row and column dimensions do not match");
}


#TODO IMPLEMENT
#SDiagMatrix is a subdiagonal matrix (zeros on and above the diagonal)
#

#########################################################################/**
# @RdocMethod setValuesWithList
#
# @title "Set Matrix Values With A List"
# 
# \description{
#   This abstract method is overwritten by all subclasses
#   to insert a list of floating-point numbers into the values
#   of a matrix.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxMatrix object.}
#  \item{valuesList}{The vector of values to insert.}
#  \item{...}{unused}
# }
#
#*/######################################################################### 
setMethodS3("setValuesWithList", "MxMatrix", 
    function(this, valuesList,...) {})


#########################################################################/**
# @RdocMethod setParametersWithList
#
# @title "Set Matrix Parameters With A List"
# 
# \description{
#   This abstract method is overwritten by all subclasses
#   to insert a list of numbers into the parameters
#   of a matrix.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxMatrix object.}
#  \item{parametersList}{The vector of parameters to insert.}
#  \item{...}{unused}
# }
#
#*/######################################################################### 
setMethodS3("setParametersWithList", "MxMatrix", 
    function(this, parametersList,...) {})


#########################################################################/**
# @RdocMethod checkValidMatrix
#
# @title "Inspect a Data Matrix for Validity"
# 
# \description{
#   This abstract method is overwritten by all subclasses
#   to return TRUE or FALSE based on whether a candidate
#   matrix can be used as a data matrix for this type.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxMatrix object.}
#  \item{aMatrix}{The candidate data matrix.}
#  \item{...}{unused}
# }
#
#*/######################################################################### 
setMethodS3("checkValidMatrix", "MxMatrix", 
    function(this, aMatrix,...) { FALSE })

#########################################################################/**
# @RdocMethod checkValidSpecification
#
# @title "Inspect a Parameter Matrix for Validity"
# 
# \description{
#   This abstract method is overwritten by all subclasses
#   to return TRUE or FALSE based on whether a candidate
#   matrix can be used as a parameter matrix for this type.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxMatrix object.}
#  \item{aMatrix}{The candidate parameter matrix.}
#  \item{...}{unused}
# }
#
#*/######################################################################### 
setMethodS3("checkValidSpecification", "MxMatrix", 
    function(this, aMatrix,...) { FALSE })


    
###########################################################################/**
# @RdocMethod updateMatrix
#
# @title "Update A Matrix Value"
# \description{
#   This is a helper function to updateMatricesHelper.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxMatrix object.}
#  \item{row}{The row of the values matrix.}
#  \item{col}{The column of the values matrix.}
#  \item{newvalue}{The new value to insert.}
#  \item{...}{unused}
# }
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("updateMatrix", "MxMatrix", function(this, row, col, newvalue, ...) {
   this$.values[row,col] <- newvalue;
});



setMethodS3("getParameters", "MxMatrix", function(this, ...) {
   this$.parameters;
});

setMethodS3("getValues", "MxMatrix", function(this, ...) {
   this$.values;
});     