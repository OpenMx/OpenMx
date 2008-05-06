#
# MxMatrix is the abstract supertype of all matrix types
#
# parameters is a matrix of integers: 
#	0    indicates a fixed value,
#	> 0  indicates a free parameter,
#       < 0  indicates FOO.
# 
# data is a matrix of floating-point values
#
# modifiable is an integer that returns the number
# of modifiable elements of this matrix type.
#
setConstructorS3("MxMatrix", function(parameters, data, modifiable) {

  if (missing(parameters)) parameters <- NA;
  if (missing(data))  data <- NA;
  if (missing(modifiable)) modifiable <- NA;

  extend(Object(), "MxMatrix",
    .parameters=parameters,
    .data=data,
    .modifiable=modifiable
  );

})


setMethodS3("print", "MxMatrix", function(this, ...) {
   cat(paste("MxMatrix:", data.class(this)), sep="\n")
   cat("Parameters: ", sep="\n")
   print(this$.parameters)
   cat("Data: ", sep="\n")
   print(this$.data)
   cat("Modifiable:", this$.modifiable)
   cat("\n")
   invisible(this)
})


checkMatrix <- function(row,col) {
   if (row < 0) throw("Row value is a negative number");
   if (col < 0) throw("Col value is a negative number");
   if (trunc(row) != row) throw("Row value is not an integer");
   if (trunc(col) != col) throw("Col value is not an integer");
}

checkFree <- function(free) {
   if (free == TRUE) warning ("Ignoring free elements on a fixed matrix");
}

checkSquare <- function(row, col) {
   if (row != col) throw("Row and column dimensions do not match");
}

#
# ZeroMatrix is full of zeros
#
setConstructorS3("ZeroMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkFree(free);

   freeParameters <- matrix(0, row, col);
   dataMatrix <- matrix(0, row, col);
   modifiable <- 0;
   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "ZeroMatrix");

})

#
# UnitMatrix is full of ones
#
setConstructorS3("UnitMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkFree(free);
   
   freeParameters <- matrix(0, row, col);
   dataMatrix <- matrix(1, row, col);
   modifiable <- 0;
   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "UnitMatrix");

})

#
# IdenMatrix is the square identity matrix
#
setConstructorS3("IdenMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkFree(free);
   checkSquare(row,col);

   freeParameters <- matrix(0, row, col);
   dataMatrix <- diag(row);
   modifiable <- 0;
   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "IdenMatrix");   

})



#
# IZeroMatrix is an Identity|Zero Partitioned Matrix
#
setConstructorS3("IZeroMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkFree(free);   

   freeParameters <- matrix(0, row, col);
   if (row == col) {
      dataMatrix <- diag(row);
   } else if (row < col) {
      square <- diag(row);
      zeros  <- matrix(0, row, col - row);
      dataMatrix <- cbind(square, zeros);
   } else {
      square <- diag(col);
      zeros  <- matrix(0, col, row - col);	
      dataMatrix <- rbind(square, zeros);
   }
   modifiable <- 0;
   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "IZeroMatrix");

})

#
# ZIdenMatrix is an Zero Partitioned Matrix|Identity
#
setConstructorS3("ZIdenMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkFree(free);

   freeParameters <- matrix(0, row, col);
   if (row == col) {
      dataMatrix <- diag(row);
   } else if (row < col) {
      square <- diag(row);
      zeros  <- matrix(0, row, col - row);
      dataMatrix <- cbind(zeros, square);
   } else {
      square <- diag(col);
      zeros  <- matrix(0, col, row - col);
      dataMatrix <- rbind(zeros, square);
   }
   modifiable <- 0;
   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "ZIdenMatrix");

})


#
# DiagMatrix is a diagonal matrix
#
setConstructorS3("DiagMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkSquare(row,col);  
  
   freeParameters <- matrix(0, row, col);
   if (free) freeParameters <- diag(row) * 1:row;
   dataMatrix <- matrix(0, row, col);
   modifiable <- row;

   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "DiagMatrix");

})

setMethodS3("checkValidMatrix", "DiagMatrix", function(this, aMatrix,...) {
   ltriangle <- all(aMatrix[lower.tri(aMatrix, diag=FALSE)] == 0);
   utriangle <- all(aMatrix[upper.tri(aMatrix, diag=FALSE)] == 0);
   return(ltriangle && utriangle);
})

setMethodS3("checkValidSpecification", "DiagMatrix", function(this, aMatrix,...) {
   ltriangle <- all(aMatrix[lower.tri(aMatrix, diag=FALSE)] == 0);
   utriangle <- all(aMatrix[upper.tri(aMatrix, diag=FALSE)] == 0);
   return(ltriangle && utriangle);
})


setMethodS3("setValuesWithList", "DiagMatrix", function(this, valuesList,...) {
   this$.data <- diag(valuesList);
})


setMethodS3("setParametersWithList", "DiagMatrix", function(this, parametersList,...) {
   this$.parameters <- diag(parametersList);
})


#
#SDiagMatrix is a subdiagonal matrix (zeros on and above the diagonal)
#



#
# StandMatrix is a standardized matrix (symmetric, ones on the diagonal)
#

setConstructorS3("StandMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkSquare(row,col);   

   dataMatrix <- diag(row);
   triangle <- lower.tri(dataMatrix, diag = FALSE);
   modifiable <- length(triangle[triangle]);
   freeParameters <- matrix(0, row, col);
   if (free) {
      ones <- matrix(1, row, col);
      # Set the lower triangular matrix to zero
      ones[lower.tri(ones, diag=TRUE)] <- 0;
      # And set the remaining elements to 1, 2, ... , n
      ones[!lower.tri(ones, diag=TRUE)] <- 1:sum(ones);
      freeParameters <- ones + t(ones);
   }

   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "StandMatrix");

})


setMethodS3("checkValidMatrix", "StandMatrix", function(this, aMatrix,...) {
   symmetry <- all(aMatrix == t(aMatrix));
   ones <- all(diag(aMatrix) == 1);
   return(symmetry && ones);
})

setMethodS3("checkValidSpecification", "StandMatrix", function(this, aMatrix,...) {
   symmetry <- all(aMatrix == t(aMatrix));
   ones <- all(diag(aMatrix) == 0);
   return(symmetry && ones);
})

setMethodS3("setValuesWithList", "StandMatrix", function(this, valuesList,...) {

   # Set the lower triangular matrix to zero
   this$.data[lower.tri(this$.data, diag=TRUE)] <- 0;

   # And set the remaining elements to valuesList
   this$.data[!lower.tri(this$.data, diag=TRUE)] <- valuesList;
   this$.data <- matrix(this$.data, nrow(this$.data), ncol(this$.data), byrow=TRUE);
   this$.data <- this$.data + t(this$.data) + diag(nrow = nrow(this$.data));
})

setMethodS3("setParametersWithList", "StandMatrix", function(this, parametersList,...) {

   # Set the lower triangular matrix to zero
   this$.parameters[lower.tri(this$.parameters, diag=TRUE)] <- 0;

   # And set the remaining elements to valuesList
   this$.parameters[!lower.tri(this$.parameters, diag=TRUE)] <- parametersList;	
   this$.parameters <- matrix(this$.parameters, nrow(this$.parameters),
					ncol(this$.parameters), byrow=TRUE);
   this$.parameters <- this$.parameters + t(this$.parameters);
})



#
# SymmMatrix is a symmetric matrix
#
setConstructorS3("SymmMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkSquare(row,col);

   dataMatrix <- matrix(0, row, col);
   triangle <- lower.tri(dataMatrix, diag = TRUE);
   modifiable <- length(triangle[triangle]);
   freeParameters <- matrix(0, row, col);
   if (free) {
      ones <- matrix(1, row, col);
      # Set the upper triangular matrix to zero
      ones[lower.tri(ones, diag=TRUE)] <- 0;
      # And set the remaining elements to 1, 2, ... , n
      ones[!lower.tri(ones, diag=TRUE)] <- 1:sum(ones);
      # And now fill in the diagonals
      biggest <- max(ones);
      diagonals <- diag(row) * (biggest + 1) : (biggest + row);
      freeParameters <- ones + t(ones) + diagonals;
   }

   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "SymmMatrix");

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
   this$.data[lower.tri(this$.data, diag=TRUE)] <- valuesList;	
   this$.data[upper.tri(this$.data, diag=FALSE)] <- this$.data[lower.tri(this$.data, diag=FALSE)];
})

setMethodS3("setParametersWithList", "SymmMatrix", function(this, parametersList,...) {

   # Set the lower triangle to parametersList
   this$.parameters[lower.tri(this$.parameters, diag=TRUE)] <- parametersList;
   this$.parameters[upper.tri(this$.parameters, diag=FALSE)] <- this$.parameters[lower.tri(this$.parameters, diag=FALSE)];
})

#
# LowerMatrix is a lower triangular matrix
#
setConstructorS3("LowerMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);
   checkSquare(row,col);

   dataMatrix <- matrix(0, row, col);
   triangle <- lower.tri(dataMatrix, diag = TRUE);
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

   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "LowerMatrix");
   
})


setMethodS3("checkValidMatrix", "LowerMatrix", function(this, aMatrix,...) {
   utriangle <- all(aMatrix[upper.tri(aMatrix, diag=FALSE)] == 0);
   return(utriangle);
})

setMethodS3("checkValidSpecification", "LowerMatrix", function(this, aMatrix,...) {
   utriangle <- all(aMatrix[upper.tri(aMatrix, diag=FALSE)] == 0);
   return(utriangle);
})

setMethodS3("setValuesWithList", "LowerMatrix", function(this, valuesList,...) {
   # Set the upper triangle to valuesList
   this$.data[upper.tri(this$.data, diag=TRUE)] <- valuesList;
   this$.data <- t(this$.data);
   this$.data[upper.tri(this$.data, diag=FALSE)] <- 0;
})

setMethodS3("setParametersWithList", "LowerMatrix", function(this, parametersList,...) {
   # Set the upper triangle to parametersList
   this$.parameters[upper.tri(this$.parameters, diag=TRUE)] <- parametersList;	
   this$.parameters <- t(this$.parameters);
   this$.parameters[upper.tri(this$.parameters, diag=FALSE)] <- 0;
})

#
# FullMatrix is a full matrix
#
setConstructorS3("FullMatrix", function(row, col, free = FALSE) {

   if (missing(row)) row <- 0;
   if (missing(col)) col <- 0;

   checkMatrix(row,col);

   dataMatrix <- matrix(0, row, col);
   modifiable <- (row * col);
   freeParameters <- matrix(0, row, col);
   if (free) freeParameters <- matrix(1 : (row * col), row, col, byrow=TRUE);

   extend(MxMatrix(freeParameters, dataMatrix, modifiable), "FullMatrix");

})

setMethodS3("setValuesWithList", "FullMatrix", function(this, valuesList,...) {
   this$.data <- t(matrix(valuesList, ncol(this$.data), nrow(this$.data)));
})

setMethodS3("setParametersWithList", "FullMatrix", function(this, parametersList,...) {
   this$.parameters <- t(matrix(parametersList, ncol(this$.parameters), nrow(this$.parameters)));
})

setMethodS3("checkValidMatrix", "FullMatrix", function(this, aMatrix,...) { TRUE })

setMethodS3("checkValidSpecification", "FullMatrix", function(this, aMatrix,...) { TRUE })
