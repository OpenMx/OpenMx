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
  \item{model}{This argument is either an MxModel object or a string. If 'model' is an MxModel object, then all elements of that model are placed in the resulting MxModel object. If 'model' is a string, then a new model is created with the string as its name.}
 	\item{...}{An arbitrary number of named entities, data sources, or MxPath objects. These will all be added or removed from the model as specified in the 'model' argument, based on the 'remove' argument.}
 	\item{manifestVars}{A list of manifest variables to be included in the model.}
 	\item{latentVars}{A list of latent variables to be included in the model.}
 	\item{remove}{logical. If TRUE, elements listed in this statement are added to the original model. If FALSE, elements listed in this statement are removed from the original model.}
 	\item{independent}{logical. If TRUE then the model is independent.}
 	\item{name}{An optional character string indicating the name of the object.}
}
 	
\details{
  The mxModel function is used to create new \link{MxModel} objects. 

  All mxModel functions must begin with a 'model' argument. If a character string is given, then the \link{MxModel} object created contains only the named entities listed in the function. If an \link{MxModel} object is given for the 'model' argument, then all elements of that object are placed in the output
  
  The 'model' argument need not be explicitly referenced, but must be used before all other arguments in this function.

  Discussion of the ellipsis. Don't forget that models may contain other models.

  All of the arguments after the ellipsis ("...") must be named arguments (i.e. foo=bar).
 	   Otherwise they will be interpreted as elements of the ellipsis list.
 	   
 	Manifest and latent variable lists.
 	
 	Independence
 	
 	Relation to other functions, specifically mxJobRun
}
 	
\value{
        Returns a new \link{MxModel} object.
}
 	
\examples{
  # Create an empy model, and place it in an object.
  model <- mxModel()
 	
  # Create a model named 'first' with one matrix
  model <- mxModel('first', mxMatrix('Full', nrow = 2,
                    ncol = 3, name = "A"))
 	
  # Add another matrix to model 'first', and rename that model 'second'
  model <- mxModel(model, mxMatrix('Iden', nrow = 3, name = "I")
                   name= "second")
 	
 	}