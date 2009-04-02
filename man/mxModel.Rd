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
  \item{model}{This argument is either an \link{MxModel} object or a string. If 'model' is an MxModel object, then all elements of that model are placed in the resulting MxModel object. If 'model' is a string, then a new model is created with the string as its name. If 'model' is not specified, a }
 	\item{...}{An arbitrary number of named entities, data sources, or MxPath objects. These will all be added or removed from the model as specified in the 'model' argument, based on the 'remove' argument.}
 	\item{manifestVars}{A list of manifest variables to be included in the model.}
 	\item{latentVars}{A list of latent variables to be included in the model.}
 	\item{remove}{logical. If TRUE, elements listed in this statement are added to the original model. If FALSE, elements listed in this statement are removed from the original model.}
 	\item{independent}{logical. If TRUE then the model is independent.}
 	\item{name}{An optional character string indicating the name of the object.}
}
 	
\details{
The mxModel function is used to create \link{MxModel} objects. Objects created by this function may be new, or may be modified versions of existing \link{MxModel} objects. To create a new \link{MxModel} object, omit or specify a character string in the 'model' argument. To create a modified version of an existing \link{MxModel} object, include this model in the 'model' argument.

Other \link{named-entities} may be added as arguments to the mxModel function, which are then added to or removed from the model specified in the 'model' argument. \link{MxAlgebra} objects, \link{MxBounds} objects, \link{MxConstraint} objects, \link{MxData} objects, \link{MxMatrix} objects, \link{MxModel} objects and objective functions may all be added in this way. \link{MxModel} objects that are included as arguments will be considered sub-models of the output model, and may be estimated separately or jointly depending on shared parameters and the 'independent' flag discussed below. Only one \link{MxData} object and one objective function may be included per model, but there are no restrictions on the number of other \link{named-entities} included in an mxModel statement.

All other arguments must be named arguments (i.e. latentVars=names), else they will be interpreted as elements of the ellipsis list. The 'manifestVars' and 'latentVars' arguments specify the names of the manifest and latent variables, respectively, for use with the \link{mxPath} function. The 'remove' argument may be used when mxModel is used to create a modified version of an existing \link{MxMatrix} object. When 'remove' is set to TRUE, the listed objects are removed from the model specified in the 'model' argument. When 'remove' is set to FALSE, the listed objects are added to the model specified in the 'model' argument. 
 	
Model independence may be specified with the 'independent' argument. If a model is independent (independent=TRUE), then the parameters of this model are not shared with any other model. An independent model may be estimated with no dependency on any other model. If a model is not independent (independent=FALSE), then this model shares parameters with one or more other models such that these models must be jointly estimated. These dependent models must be entered as arguments in another model, so that they are simultaneously optimized.
 	
To be estimated, \link{MxModel} objects must include objective functions as arguments
%(\link{mxAlgebraObjective}, \link{mxFIMLObjective}, \link{mxMLObjective} or \link{mxRAMObjective})
and executed using the \link{mxRun} function. When \link{MxData} objects are included in models, the 'type' argument of these objects may require or exclude certain objective functions, or set an objective function as default. 

\link[=Named-entity]{Named entities} in \link{MxModel} objects may be viewed and referenced by name using double brackets (model[["matrixname"]]). Slots may be referenced with the @ symbol (model@data). See the documentation for \link{Extract} and the examples in this document for more information.
}
 	
\value{
Returns a new \link{MxModel} object. \link{MxModel} objects must include an objective function to be used as arguments in \link{mxRun} functions.
}

\references{
The Mx manual can be found at http://www.vcu.edu/mx/documentation.html.
}

\seealso{
\link{MxModel} for the S4 class created by mxMatrix. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
	
\examples{
# Create an empy model, and place it in an object.
model <- mxModel()
 	
# Create a model named 'firstdraft' with one matrix
model <- mxModel('firstdraft', 
                 mxMatrix('Full', nrow = 3, ncol = 3, name = "A"))
 	
# Add other matrices to model 'firstdraft', and rename that model 'finaldraft'
model <- mxModel(model,
                 mxMatrix('Symm', nrow = 3, ncol = 3, name = "S"),
                 mxMatrix('Iden', nrow = 3, name = "F"),
                 name= "finaldraft")

# Add data to the model from an existing data frame in object 'data'
model <- mxModel(model, mxData(data, type=raw))

#View the matrix named "A" in MxModel object 'model'
model[["A"]]

#View the data associated with MxModel object 'model'
model@data

}