\name{mxAlgebra}
\alias{mxAlgebra}

\title{Create MxAlgebra Object}

\description{
   This function creates a new \link{MxAlgebra} object.
}

\arguments{
   \item{expression}{An R expression of matrix operators and matrix functions.}
   \item{name}{A character string indicating the name of the object. Optional.}
}

\value{
	Returns a new \link{MxAlgebra} object.
}

\examples{

   # the algebra B evaluates to the matrix A
   A <- mxMatrix("Full", nrow = 3, ncol = 3, name = "A")
   B <- mxAlgebra(A, name = "B")

   # the algebra C is A + B
   C <- mxAlgebra(A + B, name = "C")

   # the algebra D is sin(C)
   D <- mxAlgebra(sin(C), name = "D")
}
