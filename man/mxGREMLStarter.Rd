\name{mxGREMLStarter}
\alias{mxGREMLStarter}

\title{
Helper Function for Initializing a GREML Model
}
\description{
This function takes a dataframe or matrix and uses it to setup the 'y' and 'X' matrices for a GREML analysis.  It trims out \code{NA}s from 'X' and 'y,' and creates an \code{\link{MxFitFunctionGREML}} object that will trim out the corresponding rows and columns from 'V' at runtime.  It returns an \code{\link{MxModel}} object containing the two \code{\link[=MxMatrix]{MxMatrices}} objects and the \code{\link{MxFitFunctionGREML}} object.
}
\usage{
mxGREMLStarter(model, data, Xdata, ydata, Xname="X", yname="y", 
  addOnes=TRUE, blockByPheno=TRUE, staggerZeroes=TRUE, dropNAfromV=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Can either be a character string or an object of class \code{\link{MxModel}}.  If a character string, \code{mxGREMLStarter()} returns a new \code{\link{MxModel}} the name of which is the character string.  If an \code{\link{MxModel}}, \code{mxGREMLStarter()} returns that model with the objects the function creates added to it.}
  \item{data}{Either a dataframe or matrix, with column names, containing the variables to be used as phenotypes and covariates in 'y' and 'X,' respectively.}
  \item{Xdata}{A list of character vectors.  The number of elements (vectors) in the list equals the number of phenotypes.  The elements of each vector (the character strings) are the column names of \code{data} to be used as covariates with the corresponding phenotype.}
  \item{ydata}{Character vector.  Each string is a column name of argument \code{data}, to be used as a phenotype.}
  \item{Xname}{Character string.  The name for the \code{\link{MxMatrix}} that \code{mxGREMLStarter()} will create, to serve as 'X' in the GREML analysis.}
  \item{yname}{Character string.  The name for the \code{\link{MxMatrix}} that \code{mxGREMLStarter()} will create, to serve as 'y' in the GREML analysis.}
  \item{addOnes}{Logical; should lead columns of ones (for the regression intercepts) be adhered to the covariates when assembling the 'X' matrix?  Defaults to \code{TRUE}.}
  \item{blockByPheno}{Logical; relevant to polyphenotype analyses.  If \code{TRUE} (default), then the resulting 'y' will contain phenotype #1 for individuals 1 thru \emph{n}, phenotype #2 for individuals 1 thru \emph{n}, ...  If \code{FALSE}, then observations are "blocked by individual", and the resulting 'y' will contain individual #1's scores on phenotypes 1 thru \emph{p}, individual #2's scores on phenotypes 1 thru \emph{p}, ...  Note that 'X' is always structured appropriately for 'y.'}
  \item{staggerZeroes}{Logical; relevant to polyphenotype analyses.  If \code{TRUE} (default), then each phenotype's covariates in 'X' are "staggered," and 'X' is padded out with zeroes.  If \code{FALSE}, then 'X' is formed simply by stacking the phenotypes' covariates; this requires each phenotype to have the same number of covariates (i.e., each character vector in \code{Xdata} must be of the same length).  The default (\code{TRUE}) is intended for instances where the multiple phenotypes truly are different variables, whereas \code{staggerZeroes=FALSE} is intended for instances where the multiple "phenotypes" actually represent multiple observations on the same variable.  One example of the latter case is longitudinal data, where the multiple "phenotypes" would be repeated measures on a single phenotype.}
  \item{dropNAfromV}{Logical; internally passed to \code{\link{mxFitFunctionGREML}()}.}
}
\details{
In some cases, it may be simple for the user to drop rows and columns from 'V' prior to runtime.  In such cases, use \code{dropNAfromV=FALSE}.  Then, look at the \code{casesToDrop} slot of the fitfunction in the \code{\link{MxModel}} output by \code{mxGREMLStarter()}.  There, the indices of the rows and columns that need to be dropped can be found.  Since internally dropping rows and columns during runtime has a performance cost, it is recommended that it be done prior to runtime if possible.

For a monophenotype analysis (only), argument \code{Xdata} can be a character vector.  In a polyphenotype analysis, if the same covariates are to be used with all phenotypes, then \code{Xdata} can be a list of length 1.
}
\value{
The function returns an \code{\link{MxModel}}, which contains the \code{\link{MxFitFunctionGREML}} and the two \code{\link[=MxMatrix]{MxMatrices}} to serve as 'y' and 'X.'  As of this writing, GREML analyses do not make use of \code{\link{MxData}} objects.  But, for compatibility reasons, the returned \code{\link{MxModel}} also contains a "dummy" \code{\link{MxData}} object.
}
\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
For more information about the GREML fitfunction, see \code{\link{mxFitFunctionGREML}()} and \code{\link{MxFitFunctionGREML}}.  For more information generally concerning GREML analyses, including a complete example, see \code{\link{mxExpectationGREML}()}.  More information about the OpenMx package may be found \link[=OpenMx]{here}.
}
\examples{
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")
dat[42,1] <- NA
dat[57,2] <- NA
start <- mxGREMLStarter("GREMLmod", data=dat, Xdata="x", ydata = "y",
	addOnes = FALSE, dropNAfromV = TRUE)
str(start)
}
