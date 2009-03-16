\name{MxModel-class}
\alias{MxModel}

\title{MxModel Class}

\description{
   MxModel is an S4 class. An MxModel object is a \link[=Named-entity]{named entity}.
   New instances of this class can be created using the function \link{mxModel}.
}

\details{
  The MxModel class has the following slots:

  \tabular{rcl}{
	\tab \tab \cr
	name \tab - \tab The name of the object \cr
	matrices \tab - \tab A list of \link{MxMatrix} objects \cr
	algebras \tab - \tab A list of \link{MxAlgebra} objects \cr
	submodels \tab - \tab A list of MxModel objects \cr
	constraints \tab - \tab A list of MxConstraint objects \cr
	bounds \tab - \tab A list of MxBounds objects \cr
	paths \tab - \tab A data frame of paths \cr
	latentVars \tab - \tab A list of latent variables \cr
	manifestVars \tab - \tab A list of manifest variables \cr
	data \tab - \tab A \link{MxData} object \cr
	objective \tab - \tab Either NULL or a MxObjective object \cr
	independent \tab - \tab TRUE if-and-only-if the model is independent \cr
	output \tab - \tab A list with optimization results
  }
}
