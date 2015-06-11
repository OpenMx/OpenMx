%
%   Copyright 2007-2015 The OpenMx Project
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

\name{MxMatrix-class}
\Rdversion{1.1}
\docType{class}
\alias{MxMatrix-class}
\alias{MxMatrix}
\alias{$,MxMatrix-method}
\alias{$<-,MxMatrix-method}
\alias{[,MxMatrix-method}
\alias{[<-,MxMatrix-method}
\alias{[[,MxMatrix-method}
\alias{[[<-,MxMatrix-method}
\alias{dim,MxMatrix-method}
\alias{dimnames,MxMatrix-method}
\alias{dimnames<-,MxMatrix-method}
%\alias{imxCreateMatrix,MxMatrix-method}
%\alias{imxDeparse,MxMatrix-method}
%\alias{imxSquareMatrix,MxMatrix-method}
%\alias{imxSymmetricMatrix,MxMatrix-method}
%\alias{imxVerifyMatrix,MxMatrix-method}
\alias{DiagMatrix-class}
\alias{FullMatrix-class}
\alias{IdenMatrix-class}
\alias{LowerMatrix-class}
\alias{SdiagMatrix-class}
\alias{StandMatrix-class}
\alias{SymmMatrix-class}
\alias{UnitMatrix-class}
\alias{ZeroMatrix-class}
\alias{length,MxMatrix-method}
\alias{names,MxMatrix-method}
\alias{ncol,MxMatrix-method}
\alias{nrow,MxMatrix-method}
\alias{print,MxMatrix-method}
\alias{show,MxMatrix-method}

\title{MxMatrix Class}
\description{
MxMatrix is a virtual S4 class that comprises the nine types of matrix objects used by OpenMx (see \code{\link{mxMatrix}()} for details).  An MxMatrix object is a \link[=Named-entity]{named entity}.  New instances of this class can be created using the function \code{\link{mxMatrix}()}.  MxMatrix objects may be used as arguments in other functions from the OpenMx package, including \code{\link{mxAlgebra}()}, \code{\link{mxConstraint}()}, and \code{\link{mxModel}()}.
}
\section{Objects from the Class}{All nine types of object that the class comprises can be created via \code{\link{mxMatrix}()}.}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Character string; the name of the MxMatrix object.  Note that this is the object's "Mx name" (so to speak), which identifies it in OpenMx's internal namespace, rather than the symbol identifying it in \R's worskpace.  Use of MxMatrix objects in an \link{mxAlgebra} or \link{mxConstraint} function requires reference by name.}
    \item{\code{values}:}{Numeric matrix of values.  If an element is specified as a fixed parameter in the 'free' matrix, then the element in the 'values' matrix is treated as a constant value and cannot be altered or updated by an objective function when included in an \code{\link{mxRun}()} function. If an element is specified as a free parameter in the 'free' matrix, the element in the 'value' matrix is considered a starting value and can be changed by an objective function when included in an \code{\link{mxRun}()} function.}
    \item{\code{labels}:}{Matrix of character strings which provides the labels of free and fixed parameters. Fixed parameters with identical labels must have identical values.  Free parameters with identical labels impose an equality constraint. The same label cannot be applied to a free parameter and a fixed parameter. A free parameter with the label 'NA' implies a unique free parameter, that cannot be constrained to equal any other free parameter.}
    \item{\code{free}:}{Logical matrix specifying whether each element is free versus fixed. An element is a free parameter if-and-only-if the corresponding value in the 'free' matrix is 'TRUE'. Free parameters are elements of an MxMatrix object whose values may be changed by a fitfunction when that MxMatrix object is included in an \link{MxModel} object and evaluated using the \code{\link{mxRun}()} function.}
    \item{\code{lbound}:}{Numeric matrix of lower bounds on free parameters.}
    \item{\code{ubound}:}{Numeric matrix of upper bounds on free parameters.}
    \item{\code{.squareBrackets}:}{Logical matrix; used internally by OpenMx.  Identifies which elements have labels with square brackets in them.}
    \item{\code{.persist}:}{Logical; used internally by OpenMx.  Governs how \code{\link{mxRun}()} handles the MxMatrix object when it is inside the \link{MxModel} being run.}
    \item{\code{.condenseSlots}:}{Logical; used internally by OpenMx.  If \code{FALSE}, then the matrices in the 'values', 'labels', 'free', 'lbound', and 'ubound' slots are all of equal dimensions.  If \code{TRUE}, then the last four of those slots will "condense" a matrix consisting entirely of \code{FALSE} or \code{NA} down to 1x1.}
    \item{\code{display}:}{Character string; used internally by OpenMx when parsing \link{MxAlgebra}s.}
    \item{\code{dependencies}:}{Integer; used internally by OpenMx when parsing \link{MxAlgebra}s.}
  }
}

%Is this section required for CRAN compliance?:
\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "MxMatrix")}: ... }
    \item{$<-}{\code{signature(x = "MxMatrix")}: ... }
    \item{[}{\code{signature(x = "MxMatrix")}: ... }
    \item{[<-}{\code{signature(x = "MxMatrix")}: ... }
    \item{dim}{\code{signature(x = "MxMatrix")}: ... }
    \item{dimnames}{\code{signature(x = "MxMatrix")}: ... }
    \item{dimnames<-}{\code{signature(x = "MxMatrix")}: ... }
    %\item{imxCreateMatrix}{\code{signature(.Object = "MxMatrix")}: ... }
    %\item{imxDeparse}{\code{signature(object = "MxMatrix")}: ... }
    %\item{imxSquareMatrix}{\code{signature(.Object = "MxMatrix")}: ... }
    %\item{imxSymmetricMatrix}{\code{signature(.Object = "MxMatrix")}: ... }
    %\item{imxVerifyMatrix}{\code{signature(.Object = "MxMatrix")}: ... }
    \item{length}{\code{signature(x = "MxMatrix")}: ... }
    \item{names}{\code{signature(x = "MxMatrix")}: ... }
    \item{ncol}{\code{signature(x = "MxMatrix")}: ... }
    \item{nrow}{\code{signature(x = "MxMatrix")}: ... }
    \item{print}{\code{signature(x = "MxMatrix")}: ... }
    \item{show}{\code{signature(object = "MxMatrix")}: ... }
   }
Note that some methods are documented separately (see below, under "See Also").
}
\references{
The OpenMx User's guide can be found at \url{http://openmx.psyc.virginia.edu/documentation}.
}
\seealso{
\code{\link{mxMatrix}()} for creating MxMatrix objects.  Note that functions \code{\link{imxCreateMatrix}()}, \code{\link{imxDeparse}()}, \code{\link{imxSquareMatrix}()}, \code{\link{imxSymmetricMatrix}()}, and \code{\link{imxVerifyMatrix}()} are separately documented methods for this class.  More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}
\examples{
showClass("MxMatrix")
}
\keyword{classes}
