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
	Every named entity is guarunteed to have a slot 'name'. Within a model,
	the named entities of that model can be accessed using the model[['name']]
	notation. Access is limited to one nesting depth, such that if 'bar' is a submodel
	of 'foo', and 'baz' is a matrix of 'bar', then 'baz' must be accessed using
	foo[['bar']][['baz']].

	The following S4 classes are named entities in the OpenMx library:
	\link{MxAlgebra}, \link{MxBounds}, \link{MxConstraint}, 
	\link{MxMatrix}, \link{MxModel}, \link{MxData}, and MxObjective.

}
