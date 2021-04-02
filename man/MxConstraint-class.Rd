%
%   Copyright 2007-2018 by the individuals mentioned in the source code history
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

\name{MxConstraint-class}
\Rdversion{1.1}
\docType{class}
\alias{MxConstraint-class}
\alias{MxConstraint}
\alias{$<-,MxConstraint-method}
\alias{$,MxConstraint-method}
\alias{names,MxConstraint-method}
\alias{print,MxConstraint-method}
\alias{show,MxConstraint-method}

\title{Class \code{"MxConstraint"}}

\description{
   MxConstraint is an S4 class. An MxConstraint object
   is a \link[=Named-entity]{named entity}. New instances of this class can
   be created using the function \code{\link{mxConstraint}()}.
}

\section{Slots}{
  \describe{
    \item{\code{name}:}{Character string; the name of the object.}
    \item{\code{formula}:}{Object of class \code{"MxAlgebraFormula"}.  The \link{MxAlgebra}-like expression representing the constraint function.}
    \item{\code{alg1}:}{Object of class \code{"MxCharOrNumber"}.  For internal use.}
    \item{\code{alg2}:}{Object of class \code{"MxCharOrNumber"}.  For internal use.}
    \item{\code{relation}:}{Object of class \code{"MxCharOrNumber"}.  For internal use.}
    \item{\code{jac}:}{Object of class \code{"MxCharOrNumber"}.  Identifies the \link{MxAlgebra} representing the Jacobian for the constraint function.}
    \item{\code{linear}:}{Logical. For internal use.}
    \item{\code{strict}:}{Logical. Whether to require that all Jacobian entries reference
   free parameters.}
    \item{\code{verbose}:}{integer. For values greater than zero, enable runtime
   diagnostics.}
  }
}

\section{Methods}{
  \describe{
    \item{$<-}{\code{signature(x = "MxConstraint")}}
    \item{$}{\code{signature(x = "MxConstraint")}}
    \item{\link{imxDeparse}}{\code{signature(object = "MxConstraint")}}
    \item{names}{\code{signature(x = "MxConstraint")}}
    \item{print}{\code{signature(x = "MxConstraint")}}
    \item{show}{\code{signature(object = "MxConstraint")}}
	 }
}

\details{
Slots may be referenced with the $ symbol. See the documentation for \link[methods]{Classes} and the examples in the \link{mxConstraint} document for more information.
}

\references{
The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}
}

\seealso{
\code{\link{mxConstraint}()} for the function that creates MxConstraint objects.
}

\examples{
showClass("MxConstraint")
}

\keyword{classes}
