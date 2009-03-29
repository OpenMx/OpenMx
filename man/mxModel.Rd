\name{mxModel}
\alias{mxModel}

\title{Function To Create MxModel Object}

\description{
   This function creates a new \link{MxModel} object.
}

\usage{
mxModel(model = NA, ..., manifestVars = NA, latentVars = NA,
        remove = FALSE, independent = NA, name = NA)
}

\arguments{
   \item{model}{This argument is either an MxModel object or a string.
		If 'model' is an MxModel object, then we ether add or remove from
 		this model.  If 'model' is a string, then create a new model with
		the string as its name.}
   \item{...}{If remove = FALSE, then this argument is an arbitrary number of 
		\link[=Named-entity]{named entities} or MxPath objects. 
		If remove = TRUE, then this argument is	an arbitrary number of names
		(character strings) or MxPath objects.}
   \item{manifestVars}{A list of manifest variables}
   \item{latentVars}{A list of latent variables}
   \item{remove}{Logical. If TRUE then remove elements from the model. Defaults to FALSE.}
   \item{independent}{Logical. If TRUE then the model is independent.}
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

	# Create a model named "foo" with one matrix named "A"
	model <- mxModel("foo", mxMatrix("Full", nrow = 2, 
		ncol = 3, name = "A")) 

 	# Add matrix "I" to model "foo"
	model <- mxModel(model, mxMatrix("Iden", nrow = 3, name = "I")) 

	# Remove matrix "A" from the model
	model <- mxModel(model, "A", remove = TRUE)

}
