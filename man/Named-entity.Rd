\name{Named-entity}
\alias{Named-entity}
\alias{Named-entities}
\alias{named-entity}
\alias{named-entities}

\title{Named Entities}

\description{
	A named entity is an S4 object that can be referenced by name.
}

\details{
Every named entity is guaranteed to have a slot 'name'. Within a model, the named entities of that model can be accessed using the model[['name']] notation. Access is limited to one nesting depth, such that if 'B' is a submodel of 'A', and 'C' is a matrix of 'B', then 'C' must be accessed using A[['B']][['C']]. See the documentation for \link{Extract} for more information.

The following S4 classes are named entities in the OpenMx library: \link{MxAlgebra}, \link{MxConstraint},  \link{MxMatrix}, \link{MxModel}, \link{MxData}, and MxObjective.
}
