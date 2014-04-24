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

\name{mxThreshold}
\alias{mxThreshold}
\alias{MxThreshold-class}
\alias{print,MxThreshold-method}
\alias{show,MxThreshold-method}

\title{Create List of Thresholds}

\description{
   This function creates a list of thresholds.
}

\usage{
mxThreshold(vars, nThresh=NA,
	free=FALSE, values=NA, labels=NA,
	lbound=NA, ubound=NA)
}

\arguments{
   \item{vars}{character vector. These are the variables for which thresholds are to be specified.}
   \item{nThresh}{numeric vector. These are the number of thresholds for each variables listed in \sQuote{vars}.}
   \item{free}{boolean vector. Indicates whether threshold parameters are free or fixed.}
   \item{values}{numeric vector. The starting values of the parameters.}
   \item{labels}{character vector. The names of the parameters.}
   \item{lbound}{numeric vector. The lower bounds of free parameters.}
   \item{ubound}{numeric vector. The upper bounds of free parameters.}
}

\details{
The mxPath function creates \link[=MxThreshold-class]{MxThreshold} objects. These consist of a list of ordinal variables and the thresholds that define the relationship between the observed ordinal variable and the continuous latent variable assumed to underly it. This function directly mirrors the usage of \link{mxPath}, but is used to specify thresholds rather than means, variances and bivariate relationships.

The \sQuote{vars} argument specifies which variables you wish to specify thresholds for. Variables are referenced by name, and these names must appear in the \sQuote{manifestVar} argument of the \link{mxModel} function if thresholds are to be correctly processed. Additionally, variables for which thresholds are specified must be specified as ordinal factors in whatever data is included in the model.

The \sQuote{nThresh} argument specifies how many thresholds are to be specified for the variable or variables included in the \sQuote{vars} argument. The number of thresholds for a particular variable should be one fewer than the number of categories specified for that variable. 
    
The \sQuote{free} argument specifies whether the thresholds created by the mxThreshold function are free or fixed parameters. This argument may take either TRUE for free parameters, FALSE for fixed parameters, or a vector of TRUEs and FALSEs to be applied in order to the created thresholds.

The \sQuote{values} is a numeric vectors containing the starting values of the created thresholds. \sQuote{values} gives a starting value for estimation. The \sQuote{labels} argument specifies the names of the parameters in the resulting \link[=MxThreshold-class]{MxThreshold} object.  The \sQuote{lbound} and \sQuote{ubound} arguments specify lower and upper bounds for the created threshold parameters.

Thresholds for multiple variables may be specified simultaneously by including a vector of variable names to the \sQuote{vars} argument. When multiple variables are included in the \sQuote{vars} argument, the length of the \sQuote{vars} argument must be evenly divisable by the length of the \sQuote{nThresh} argument. All subsequent arguments (\sQuote{free} through \sQuote{ubound}) should have their lengths be a factor of the total number of thresholds specified for all variables. 

If four variables are included in the \sQuote{vars} argument, then the \sQuote{nThresh} argument should contain either one, two or four elements. If the \sQuote{nThresh} argument specifies two thresholds for each variable, then \sQuote{free}, \sQuote{values}, and all subsequent arguments should specify eight values by including one, two, four or eight elements. Whenever fewer values are specified than are required (e.g., specify two values for eight thresholds), then the entire vector of values is repeated until the required number of values is reached, and will return an error if the correct number of values cannot be achieved by repeating the entire vector.
}

\value{
Returns a list of thresholds.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxPath} for comparable specification of paths. \link{mxMatrix} for a matrix-based approach to thresholds specification; \link{mxModel} for the container in which mxThresholds are embedded. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}


\examples{
# a simple one factor ordinal model
require(OpenMx)
	
data(myFADataRaw)

oneFactorOrd <- myFADataRaw[,c("z1", "z2", "z3")]

oneFactorOrd$z1 <- mxFactor(oneFactorOrd$z1, levels=c(0, 1))
oneFactorOrd$z2 <- mxFactor(oneFactorOrd$z2, levels=c(0, 1))
oneFactorOrd$z3 <- mxFactor(oneFactorOrd$z3, levels=c(0, 1, 2))

	oneFactorModel <- mxModel("Common Factor Model Path Specification", 
	type="RAM",
	mxData(
		observed=oneFactorOrd,
		type="raw"
	),
	manifestVars=c("z1","z2","z3"),
	latentVars="F1",
	# residual variances
	mxPath(
		from=c("z1","z2","z3"),
		arrows=2,
		free=FALSE,
		values=c(1,1,1),
		labels=c("e1","e2","e3")
	),
	# latent variance
	mxPath(
		from="F1",
		arrows=2,
		free=TRUE,
		values=1,
		labels ="varF1"
	),
	# factor loadings
	mxPath(
		from="F1",
		to=c("z1","z2","z3"),
		arrows=1,
		free=c(FALSE,TRUE,TRUE),
		values=c(1,1,1),
		labels=c("l1","l2","l3")
	),
	# means
	mxPath(
		from="one",
		to=c("z1","z2","z3","F1"),
		arrows=1,
		free=FALSE,
		values=0,
		labels=c("meanz1","meanz2","meanz3","meanF")
	),
	# thresholds
	mxThreshold(vars=c("z1", "z2", "z3"),
		nThresh=c(1,1,2),
		free=TRUE,
		values=c(-1, 0, -.5, 1.2)
		)
) # close model

oneFactorResults <- mxRun(oneFactorModel)

}
