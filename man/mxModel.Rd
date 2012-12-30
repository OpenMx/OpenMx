%
%   Copyright 2007-2012 The OpenMx Project
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

\name{mxModel}
\alias{mxModel}

\title{Create MxModel Object}
 	
\description{
   This function creates a new \link{MxModel} object.
}
 	
\usage{
mxModel(model = NA, ..., manifestVars = NA, latentVars = NA,
 	      remove = FALSE, independent = NA, type = NA, name = NA)
}
 	
\arguments{
  \item{model}{This argument is either an \link{MxModel} object or a string. If 'model' is an MxModel object, then all elements of that model are placed in the resulting MxModel object. If 'model' is a string, then a new model is created with the string as its name. If 'model' is either unspecified or 'model' is a named entity, data source, or MxPath object, then a new model is created.}
 	\item{...}{An arbitrary number of \link{mxMatrix}, \link{mxPath}, \link{mxData}, and other functions such as \link{mxConstraint}s and \link{mxCI}. These will all be added or removed from the model as specified in the 'model' argument, based on the 'remove' argument.}
 	\item{manifestVars}{For RAM-type models, A list of manifest variables to be included in the model.}
 	\item{latentVars}{For RAM-type models, A list of latent variables to be included in the model.}
 	\item{remove}{logical. If TRUE, elements listed in this statement are removed from the original model. If FALSE, elements listed in this statement are added to the original model.}
 	\item{independent}{logical. If TRUE then the model is evaluated independently of other models.}
 	\item{type}{character vector. The model type to assign to this model. Defaults to options("mxDefaultType"). See below for valid types}
 	\item{name}{An optional character vector indicating the name of the object.}
}
 	
\details{
The mxModel function is used to create \link{MxModel} objects. Objects created by this function may be new, or may be modified versions of existing \link{MxModel} objects. By default a new \link{MxModel} object will be created: To create a modified version of an existing \link{MxModel} object, include this model in the 'model' argument.

Other \link{named-entities} may be added as arguments to the mxModel function, which are then added to or removed from the model specified in the \sQuote{model} argument. Other functions you can use to add objects to the model to this way are \link{mxCI}, \link{mxAlgebra}, \link{mxBounds}, \link{mxConstraint}, \link{mxData}, and \link{mxMatrix} objects, as well as objective functions. You can also include \link{MxModel} objects as sub-models of the output model, and may be estimated separately or jointly depending on shared parameters and the \sQuote{independent} flag discussed below. Only one \link{MxData} object and one objective function may be included per model, but there are no restrictions on the number of other \link{named-entities} included in an mxModel statement.

All other arguments must be named (i.e. \sQuote{latentVars = names}), or they will be interpreted as elements of the ellipsis list. The \sQuote{manifestVars} and \sQuote{latentVars} arguments specify the names of the manifest and latent variables, respectively, for use with the \link{mxPath} function. The \sQuote{remove} argument may be used when mxModel is used to create a modified version of an existing \link{MxMatrix} object. When \sQuote{remove} is set to TRUE, the listed objects are removed from the model specified in the \sQuote{model} argument. When \sQuote{remove} is set to FALSE, the listed objects are added to the model specified in the \sQuote{model} argument. 
 	
Model independence may be specified with the \sQuote{independent} argument. If a model is independent (\sQuote{independent = TRUE}), then the parameters of this model are not shared with any other model. An independent model may be estimated with no dependency on any other model. If a model is not independent (\sQuote{independent = FALSE}), then this model shares parameters with one or more other models such that these models must be jointly estimated. These dependent models must be entered as arguments in another model, so that they are simultaneously optimized.

The model type is determined by a character vector supplied to the \sQuote{type} argument. The type of a model is a dynamic property, ie. it is allowed to change during the lifetime of the model. To see a list of available types, use the \link{mxTypes} command.  When a new model is created and no type is specified, the type specified by \code{options("mxDefaultType")} is used.
 	
To be estimated, \link{MxModel} objects must include objective functions as arguments
(\link{mxAlgebraObjective}, \link{mxFIMLObjective}, \link{mxMLObjective} or \link{mxRAMObjective})
and executed using the \link{mxRun} function. When \link{MxData} objects are included in models, the 'type' argument of these objects may require or exclude certain objective functions, or set an objective function as default. 

\link[=Named-entity]{Named entities} in \link{MxModel} objects may be viewed and referenced by name using double brackets (\code{model[["matrixname"]]}). Slots may be referenced with the @ symbol (\code{model@data}). See the documentation for \link[methods]{Classes} and the examples in this document for more information.
}
 	
\value{
Returns a new \link{MxModel} object. \link{MxModel} objects must include an objective function to be used as arguments in \link{mxRun} functions.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
See \link{mxCI} for information about adding Confidence Interval calculations to a model. 
See \link{mxPath} for information about adding paths to RAM-type models.
See \link{mxMatrix} for information about adding matrices to models.
See \link{mxData} for specifying the data a model is to be evaluated against.
See \link{MxModel} for the S4 class created by mxMatrix. Many advanced options can be set via \link{mxOption} (such as calculating the Hessian). 
More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
	
\examples{
# At the simplest, you can create an empty model, placing it in an object, and add to it later
model <- mxModel()

# Create a model named 'firstdraft' with one matrix 'A'
model <- mxModel('firstdraft', 
                 mxMatrix('Full', nrow = 3, ncol = 3, name = "A"))

# Update 'firstdraft', and rename the model 'finaldraft'
model <- mxModel(model,
                 mxMatrix('Symm', nrow = 3, ncol = 3, name = "S"),
                 mxMatrix('Iden', nrow = 3, name = "F"),
                 name= "finaldraft")

# Add data to the model from an existing data frame in object 'data'
data(twinData)  # load some data
model <- mxModel(model, mxData(twinData, type='raw'))

# View the matrix named "A" in MxModel object 'model'
model[["A"]]

model@matrices$A

#View the data associated with MxModel object 'model'
model$data

# A working example using OpenMx Path Syntax
data(HS.fake.data)  #load the data

Spatial   <- c("visual","cubes","paper") # the manifest variables loading on each proposed latent variable
Verbal    <- c("general","paragrap","sentence")
Math      <- c("numeric","series","arithmet")

latents   <- c("vis","math","text")
manifests <-  c(Spatial,Math,Verbal)

model <- mxModel("Holzinger and Swineford (1939)", type="RAM", 
	manifestVars = manifests, # list the measured variables (boxes)
	latentVars   = latents,   # list the latent variables (circles)
    # factor loadings from latents to  manifests
	mxPath(from="vis",  to=Spatial),# factor loadings
    mxPath(from="math", to=Math),   # factor loadings
    mxPath(from="text", to=Verbal), # factor loadings

	# Allow latent variables to covary 
	mxPath(from="vis" , to="math", arrows=2, free=TRUE),
	mxPath(from="vis" , to="text", arrows=2, free=TRUE),
	mxPath(from="math", to="text", arrows=2, free=TRUE),

    # Allow latent variables to have variance (first fixed  @ 1)
	mxPath(from=latents, arrows=2, free=c(FALSE,TRUE,TRUE), values=1.0), 
    # Manifest have residual variance
	mxPath(from=manifests, arrows=2),   
	# the data to be analysed
    mxData(cov(HS.fake.data[,manifests]), type="cov", numObs=301))
	
fit <- mxRun(model) # run the model
summary(fit) # examine the output: Fit statistics and path loadings

}
