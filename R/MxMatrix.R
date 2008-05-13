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

#########################################################################/**
# @RdocMethod FIXED
#
# @title "Value Of A Fixed Parameter"
# 
# \description{
#   The return value is the value that is being used in all parameter
#   matrices to represent a fixed value.
# }
#
# @synopsis
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("FIXED", "MxMatrix", function(this, ...) {
   return("0");
}, static = TRUE);

#########################################################################/**
# @RdocMethod UNIQUE
#
# @title "String Prefix Of A Free Parameter"
# 
# \description{
#   The return value is the string prefix that is being used in all parameter
#   matrices to represent a unique free value.
# }
#
# @synopsis
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("UNIQUE", "MxMatrix", function(this, ...) {
   return("auto-");
}, static = TRUE);


#########################################################################/**
# @RdocMethod createUniqueName
#
# @title "Create A Unique Name For A Free Parameter"
# 
# \description{
#   The return value is garuanteed to be a unique string.
# }
#
# @synopsis
#
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
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

#########################################################################/**
# @RdocFunction checkMatrix
#
# @title "Check Matrix Dimensions"
# 
# \description{
#    A helper function that performs some sanity checking
#    on the dimensions of a new MxMatrix object.
# }
#
# @synopsis
#
# \arguments{
#  \item{row}{The proposed row size.}
#  \item{col}{The proposed column size.}
# }
#
#*/######################################################################### 
checkMatrix <- function(row,col) {
   if (row < 0) throw("Row value is a negative number");
   if (col < 0) throw("Col value is a negative number");
   if (trunc(row) != row) throw("Row value is not an integer");
   if (trunc(col) != col) throw("Col value is not an integer");
}

#########################################################################/**
# @RdocFunction checkSquare
#
# @title "Check Square Matrix Property"
# 
# \description{
#    A helper function that performs some sanity checking
#    on the dimensions of a new square MxMatrix object.
# }
#
# @synopsis
#
# \arguments{
#  \item{row}{The proposed row size.}
#  \item{col}{The proposed column size.}
# }
#
#*/######################################################################### 
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

###########################################################################/**
# @RdocMethod getParameters
#
# @title "Parameters Getter"
# \description{
#   Returns the parameters matrix which is a matrix of strings.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxMatrix object.}
#  \item{...}{unused}
# }
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getParameters", "MxMatrix", function(this, ...) {
   this$.parameters;
});

###########################################################################/**
# @RdocMethod getValues
#
# @title "Values Getter"
# \description{
#   Returns the values matrix which is a matrix of doubles.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxMatrix object.}
#  \item{...}{unused}
# }
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getValues", "MxMatrix", function(this, ...) {
   this$.values;
});     

###########################################################################/**
# @RdocMethod setValues
#
# @title "Values Setter"
# \description{
#   Sets the values matrix which is a matrix of doubles.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxMatrix object.}
#  \item{values}{The new values as a matrix or a vector.}
#  \item{...}{unused}
# }
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
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

###########################################################################/**
# @RdocMethod setParameters
#
# @title "Parameters Setter"
# \description{
#   Sets the parameter matrix which is a matrix of strings.
# }
#
# @synopsis
#
# \arguments{
#  \item{this}{The MxMatrix object.}
#  \item{parameters}{The new parameters as a matrix or a vector.}
#  \item{...}{unused}
# }
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
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
