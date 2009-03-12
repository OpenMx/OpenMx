\name{mxModel}
\alias{mxModel}

\title{Create MxModel Object}

\description{
   This function creates a new \link{MxModel} object.
}

\arguments{
   \item{model}{This argument is either an MxModel object or a string.
		If 'model' is an MxModel object, then we ether add or remove from
 		this model.  If 'model' is a string, then create a new model with
		the string as its name.}
   \item{...}{An arbitrary number of named entities, data sources, or MxPath objects.
		These will all be added or removed from the model, based on the 'remove' argument.}
   \item{manifestVars}{A list of manifest variables}
   \item{latentVars}{A list of latent variables}
   \item{remove}{Logical. If TRUE then remove elements from the model. Defaults to FALSE.}
   \item{independent}{Logical. If TRUE then the model is independent. Defaults to FALSE.}
   \item{name}{A character string indicating the name of the object. Optional.}
}

\details{
   All of the arguments after the ellipsis ("...") must be named arguments (i.e. foo=bar).
   Otherwise they will be interpreted as elements of the ellipsis list.
}

\value{
	Returns a new \link{MxModel} object.
}

\examples{

	# Create an anonymous model
	model <- mxModel()

	# Create a model named 'foo' with one matrix
	model <- mxModel('foo', mxMatrix('Full', nrow = 2, 
		ncol = 3, name = "A")) 

 	# Add another matrix to model 'foo'
	model <- mxModel(model, mxMatrix('Iden', nrow = 3, name = "I")) 

}
