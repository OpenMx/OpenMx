%
%   Copyright 2007-2019 by the individuals mentioned in the source code history
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
   Create or modify an \link{MxModel}.
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
The mxModel function is used to create \link{MxModel}s. Models created by this function may be new, or may be modified versions of existing \link{MxModel} objects. By default a new \link{MxModel} object will be created: To create a modified version of an existing \link{MxModel} object, include this model in the 'model' argument.

Other \link{named-entities} may be added as arguments to the mxModel function, which are then added to or removed from the model specified in the \sQuote{model} argument. Functions you can use to add objects to the model to this way include \link{mxPath}, \link{mxCI}, \link{mxAlgebra}, \link{mxBounds}, \link{mxConstraint}, \link{mxData}, and \link{mxMatrix} objects, as well as fit functions and expectations (see below). You can also include sub-models as components of a model. These sub-models may be estimated separately or jointly depending on shared parameters and the \sQuote{independent} flag (see below). Only one \link{MxData} object and one fit function and expectation may be included per model, but there are no restrictions on the number of other \link{named-entities} included in an mxModel statement.

All other arguments must be named (i.e. \sQuote{latentVars = names}), or they will be interpreted as elements of the ellipsis list. The \sQuote{manifestVars} and \sQuote{latentVars} arguments specify the names of the manifest and latent variables, respectively, for use with the \link{mxPath} function.

The \sQuote{remove} argument may be used when mxModel is used to create a modified version of an existing \link{MxMatrix} object. When \sQuote{remove} is set to TRUE, the listed objects are removed from the model specified in the \sQuote{model} argument. When \sQuote{remove} is set to FALSE, the listed objects are added to the model specified in the \sQuote{model} argument. 
    
Model independence may be specified with the \sQuote{independent} argument. If a model is independent (\sQuote{independent = TRUE}), then the parameters of this model are not shared with any other model. An independent model may be estimated with no dependency on any other model. If a model is not independent (\sQuote{independent = FALSE}), then this model shares parameters with one or more other models such that these models must be jointly estimated. These dependent models must be entered as arguments in another model, so that they are simultaneously optimized.

The model type is determined by a character vector supplied to the \sQuote{type} argument. The type of a model is a dynamic property, ie. it is allowed to change during the lifetime of the model. To see a list of available types, use the \link{mxTypes} command.  When a new model is created and no type is specified, the type specified by \code{options("mxDefaultType")} is used.

\strong{Expectations and Fit functions}

To be estimated, \link{MxModel} objects must include fit functions and expectations as arguments.  Fit functions include \link{mxFitFunctionML}, \link{mxFitFunctionMultigroup}, \link{mxFitFunctionWLS}, \link{mxFitFunctionAlgebra}, \link{mxFitFunctionGREML}, \link{mxFitFunctionR}, and \link{mxFitFunctionRow}. Expectations include 
\link{mxExpectationBA81}, \link{mxExpectationGREML}, \link{mxExpectationHiddenMarkov}
\link{mxExpectationLISREL}, \link{mxExpectationMixture}, \link{mxExpectationNormal}, \link{mxExpectationRAM}, \link{mxExpectationStateSpace}, \link{mxExpectationStateSpaceContinuousTime}. The 'type' of the model may imply a certain fit function or expectation (e.g. type = "RAM" implies mxExpectationRAM). The model data may also constrain which fit and expectation are appropriate.

The model, complete with fit function and expectation can then be executed using \link{mxRun}. 

\strong{Accessing model components}

You can view a model summary with summary. You can also access \link[=Named-entity]{Named entities} in \link{MxModel} directly via the $ symbol.  For instance, for an MxModel named "yourModel" containing an MxMatrix named "yourMatrix", the contents of "yourMatrix" can be accessed as yourModel$yourMatrix. Slots (i.e., matrices, algebras, etc.) in an mxMatrix may also be referenced with the $ symbol (e.g., yourModel$matrices or yourModel$algebras). See the documentation for \link[methods]{Classes} and the examples in \link[methods]{Classes} for more information.
}
    
\value{
Returns a new \link{MxModel} object. To be run, \link{MxModel} object must include a fit function and expectation.
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\seealso{
See \link{mxCI} for information about adding Confidence Interval calculations to a model. 
See \link{mxPath} for information about adding paths to RAM-type models.
See \link{mxMatrix} for information about adding matrices to models.
See \link{mxData} for specifying the data a model is to be evaluated against.
Many advanced options can be set via \link{mxOption}. 
More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
    
\examples{

library(OpenMx)

# At the simplest, you can create an empty model,
#  placing it in an object, and add to it later
emptyModel <- mxModel(model="IAmEmpty")

# Create a model named 'firstdraft' with one matrix 'A'
firstModel <- mxModel(model='firstdraft', 
                 mxMatrix(type='Full', nrow = 3, ncol = 3, name = "A"))

# Update 'firstdraft', and rename the model 'finaldraft'
finalModel <- mxModel(model=firstModel,
                 mxMatrix(type='Symm', nrow = 3, ncol = 3, name = "S"),
                 mxMatrix(type='Iden', nrow = 3, name = "F"),
                 name= "finaldraft")

# Add data to the model from an existing data frame in object 'data'
data(twinData)  # load some data
finalModel <- mxModel(model=finalModel, mxData(twinData, type='raw'))

# Two ways to view the matrix named "A" in MxModel object 'model'

finalModel$A

finalModel$matrices$A

# A working example using OpenMx Path Syntax
data(HS.ability.data)  # load the data

# The manifest variables loading on each proposed latent variable
Spatial   <- c("visual", "cubes", "paper")
Verbal    <- c("general", "paragrap", "sentence")
Math      <- c("numeric", "series", "arithmet")

latents   <- c("vis", "math", "text")
manifests <-  c(Spatial, Math, Verbal)

HSModel <- mxModel(model="Holzinger_and_Swineford_1939", type="RAM", 
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

    # Allow latent variables to have variance
    mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
    # Manifest have residual variance
    mxPath(from=manifests, arrows=2),   
    # the data to be analysed
    mxData(cov(HS.ability.data[,manifests]), type = "cov", numObs = 301))
    
fitModel <- mxRun(HSModel) # run the model
summary(fitModel) # examine the output: Fit statistics and path loadings

}
