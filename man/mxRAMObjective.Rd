\name{mxRAMObjective}
\alias{mxRAMObjective}

\title{Function To Create MxRAMObjective Object}

\description{
   This function creates a new MxRAMObjective object.
}

\usage{
mxRAMObjective(aMatrix = "A", sMatrix = "S", fMatrix = "F", name = NA)
}


\arguments{
   \item{aMatrix}{A character string indicating the name of the 'A' matrix.}
   \item{sMatrix}{A character string indicating the name of the 'S' matrix.}
   \item{fMatrix}{A character string indicating the name of the 'F' matrix.}
   \item{name}{A character string indicating the
   		name of this object. Optional.}
}

\value{
	Returns a new MxRAMObjective object.
}

\examples{
   
   A <- mxMatrix("Full", c(0,0.2,0,0), name="A", nrow=2, ncol=2)
   S <- mxMatrix("Full", c(0.8,0,0,0.8), name="S", nrow=2, ncol=2, free=TRUE)
   F <- mxMatrix("Full", c(1,0,0,1), name="F", nrow=2, ncol=2)

   # Create a RAM objective with default A, S, F matrix names
   objective <- mxRAMObjective(name = "objective")
   
   model <- mxModel(A, S, F, objective)

}