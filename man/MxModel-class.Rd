%
%   Copyright 2007-2014 The OpenMx Project
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

\name{MxModel-class}
\alias{MxModel}
\alias{$,MxModel-method}
\alias{$<-,MxModel-method}
\alias{[[,MxModel-method}
\alias{[[<-,MxModel-method}
\alias{names,MxModel-method}
\alias{print,MxModel-method}
\alias{show,MxModel-method}

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
    intervals \tab - \tab A list of confidence intervals requested in MxCI objects \cr
    bounds \tab - \tab A list of MxBounds objects \cr
    latentVars \tab - \tab A list of latent variables \cr
    manifestVars \tab - \tab A list of manifest variables \cr
    data \tab - \tab A \link{MxData} object \cr
    objective \tab - \tab Either NULL or a MxObjective object \cr
    independent \tab - \tab TRUE if-and-only-if the model is independent \cr
    options \tab - \tab A list of optimizer options \cr
    output \tab - \tab A list with optimization results
  }
  
The \sQuote{name} slot is the name of the MxModel object.

The \sQuote{matrices} slot contains a list of the \link{MxMatrix} objects included in the model. These objects are listed by name. Two objects may not share the same name. If a new \link{MxMatrix} is added to an MxModel object with the same name as an \link{MxMatrix} object in that model, the added version replaces the previous version. There is no imposed limit on the number of \link{MxMatrix} objects that may be added here.

The \sQuote{algebras} slot contains a list of the \link{MxAlgebra} objects included in the model. These objects are listed by name. Two objects may not share the same name. If a new \link{MxAlgebra} is added to an MxModel object with the same name as an \link{MxAlgebra} object in that model, the added version replaces the previous version. All \link{MxMatrix} objects referenced in the included \link{MxAlgebra} objects must be included in the \sQuote{matrices} slot prior to estimation. There is no imposed limit on the number of \link{MxAlgebra} objects that may be added here.

The \sQuote{submodels} slot contains references to all of the MxModel objects included as submodels of this MxModel object. Models held as arguments in other models are considered to be submodels. These objects are listed by name. Two objects may not share the same name. If a new submodel is added to an MxModel object with the same name as an existing submodel, the added version replaces the previous version. When a model containing other models is executed using \link{mxRun}, all included submodels are executed as well. If the submodels are dependent on one another, they are treated as one larger model for purposes of estimation.

The \sQuote{constraints} slot contains a list of the \link{MxConstraint} objects included in the model. These objects are listed by name. Two objects may not share the same name. If a new \link{MxConstraint} is added to an MxModel object with the same name as an \link{MxConstraint} object in that model, the added version replaces the previous version. All \link{MxMatrix} objects referenced in the included \link{MxConstraint} objects must be included in the \sQuote{matrices} slot prior to estimation. There is no imposed limit on the number of \link{MxAlgebra} objects that may be added here.

The \sQuote{intervals} slot contains a list of the confidence intervals requested by included \link{MxCI} objects. These objects are listed by the free parameters, \link[=MxMatrix]{MxMatrices} and \link[=MxAlgebra]{MxAlgebras} referenced in the \link{MxCI} objects, not the list of \link{MxCI} objects themselves. If a new \link{MxCI} object is added to an MxModel object referencing one or more free parameters \link[=MxMatrix]{MxMatrices} or \link[=MxAlgebra]{MxAlgebras} previously listed in the \sQuote{intervals} slot, the new confidence interval(s) replace the existing ones. All listed confidence intervals must refer to free parameters \link[=MxMatrix]{MxMatrices} or \link[=MxAlgebra]{MxAlgebras} in the model.

The \sQuote{bounds} slot contains a list of the \link{MxBounds} objects included in the model. These objects are listed by name. Two objects may not share the same name. If a new \link{MxBounds} is added to an MxModel object with the same name as an \link{MxBounds} object in that model, the added version replaces the previous version. All \link{MxMatrix} objects referenced in the included \link{MxBounds} objects must be included in the \sQuote{matrices} slot prior to estimation. There is no imposed limit on the number of \link{MxAlgebra} objects that may be added here.

The \sQuote{latentVars} slot contains a list of latent variable names, which may be referenced by \link[=MxPath-class]{MxPath} objects. This slot defaults to 'NA', and is only used when the \link{mxPath} function is used.

The \sQuote{manifestVars} slot contains a list of latent variable names, which may be referenced by \link[=MxPath-class]{MxPath} objects. This slot defaults to 'NA', and is only used when the \link{mxPath} function is used.

The \sQuote{data} slot contains an \link{MxData} object. This slot must be filled prior to execution when an objective function referencing data is used. Only one \link{MxData} object may be included per model, but submodels may have their own data in their own \sQuote{data} slots. If an \link{MxData} object is added to an MxModel which already contains an \link{MxData} object, the new object replaces the existing one.

The \sQuote{objective} slot contains an objective function. This slot must be filled prior to using the \link{mxRun} function for model execution and optimization. \link{MxAlgebra}, \link{MxData}, and \link{MxMatrix} objects required by the included objective function must be included in the appropriate slot of the MxModel prior to using \link{mxRun}.

The \sQuote{independent} slot contains a logical value indicating whether or not the model is independent. If a model is independent (independent=TRUE), then the parameters of this model are not shared with any other model. An independent model may be estimated with no dependency on any other model. If a model is not independent (independent=FALSE), then this model shares parameters with one or more other models such that these models must be jointly estimated. These dependent models must be entered as submodels of another MxModel objects, so that they are simultaneously optimized.

The \sQuote{options} slot contains a list of options for the optimizer. The name of each entry in the list is the option name to be passed to the optimizer. The values in this list are the values of the optimizer options. The standard interface for updating options is through the \link{mxOption} function.

The \sQuote{output} slot contains a list of output added to the model by the \link{mxRun} function. Output includes parameter estimates, optimization information, model fit, and other information as dictated by the objective function. If a model has not been optimized using the \link{mxRun} function, the 'output' slot will be 'NULL'.

\link[=Named-entity]{Named entities} in \link{MxModel} objects may be viewed and referenced by name using the $ symbol.  For instance, for an MxModel named "yourModel" containing an MxMatrix named "yourMatrix", the contents of "yourMatrix" can be accessed as yourModel$yourMatrix. Slots (i.e., matrices, algebras, etc.) in an mxMatrix may also be referenced with the $ symbol (e.g., yourModel$matrices or yourModel$algebras). See the documentation for \link[methods]{Classes} and the examples in \link{mxModel} for more information.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxModel} for creating MxModel objects. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
