\name{mxGREMLDataHandler}
\alias{mxGREMLDataHandler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper Function for Structuring GREML Data 
}
\description{
This function takes a dataframe or matrix and uses it to setup the 'y' and 'X' matrices for a GREML analysis; this includes trimming out \code{NA}s from 'X' and 'y.'  The result is a matrix the first column of which is the 'y' vector, and the remaining columns of which constitute 'X.'
}
\usage{
mxGREMLDataHandler(data, yvars=character(0), Xvars=list(), addOnes=TRUE, blockByPheno=TRUE, staggerZeroes=TRUE)
}

\arguments{
 \item{data}{Either a dataframe or matrix, with column names, containing the variables to be used as phenotypes and covariates in 'y' and 'X,' respectively.}
  \item{yvars}{Character vector.  Each string names a column of the raw dataset, to be used as a phenotype.}
  \item{Xvars}{A list of data column names, specifying the covariates to be used with each phenotype.  The list should have the same length as argument \code{yvars}.}
   \item{addOnes}{Logical; should lead columns of ones (for the regression intercepts) be adhered to the covariates when assembling the 'X' matrix?  Defaults to \code{TRUE}.}
   \item{blockByPheno}{Logical; relevant to polyphenotype analyses.  If \code{TRUE} (default), then the resulting 'y' will contain phenotype #1 for individuals 1 thru \emph{n}, phenotype #2 for individuals 1 thru \emph{n}, ...  If \code{FALSE}, then observations are "blocked by individual", and the resulting 'y' will contain individual #1's scores on phenotypes 1 thru \emph{p}, individual #2's scores on phenotypes 1 thru \emph{p}, ...  Note that in either case, 'X' will be structured appropriately for 'y.'}
   \item{staggerZeroes}{Logical; relevant to polyphenotype analyses.  If \code{TRUE} (default), then each phenotype's covariates in 'X' are "staggered," and 'X' is padded out with zeroes.  If \code{FALSE}, then 'X' is formed simply by stacking the phenotypes' covariates; this requires each phenotype to have the same number of covariates (i.e., each character vector in \code{Xvars} must be of the same length).  The default (\code{TRUE}) is intended for instances where the multiple phenotypes truly are different variables, whereas \code{staggerZeroes=FALSE} is intended for instances where the multiple "phenotypes" actually represent multiple observations on the same variable.  One example of the latter case is longitudinal data where the multiple "phenotypes" are repeated measures on a single phenotype.}
}
\details{
For a monophenotype analysis (only), argument \code{Xdata} can be a character vector.  In a polyphenotype analysis, if the same covariates are to be used with all phenotypes, then \code{Xdata} can be a list of length 1.

Note the synergy between the output of \code{mxGREMLDataHandler()} and arguments \code{dataset.is.yX} and \code{casesToDropFromV} to \code{\link{mxExpectationGREML}()}.

If the dataframe or matrix supplied for argument \code{data} has \emph{n} rows, and argument \code{yvars} is of length \emph{p}, then the resulting 'y' and 'X' matrices will have \emph{np} rows.  Then, if either matrix contains any \code{NA}'s, the rows containing the \code{NA}'s are trimmed from both 'X' and 'y' before being returned in the output (in which case they will obviously have fewer than \emph{np} rows).  Function \code{mxGREMLDataHandler()} reports which rows of the full-size 'X' and 'y' were trimmed out due to missing observations.  These row indices can be provided as argument \code{casesToDropFromV} to \code{\link{mxExpectationGREML}()}.

}

\value{
A list with these two components:
  \item{yX}{Numeric matrix.  The first column is the phenotype vector, 'y,' while the remaining columns constitutethe 'X' matrix of covariates.  If this matrix is used as the raw dataset for a model, then the model's GREML expectation can be constructed with \code{dataset.is.yX=TRUE} in \code{\link{mxExpectationGREML}()}.}
  \item{casesToDrop}{Numeric vector.  Contains the indices of the rows of the 'y' and 'X' that were dropped due to containing \code{NA}'s.  Can be provided as as argument \code{casesToDropFromV} to \code{\link{mxExpectationGREML}()}.}
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
For more information generally concerning GREML analyses, including a complete example, see \code{\link{mxExpectationGREML}()}.  More information about the OpenMx package may be found \link[=OpenMx]{here}.
}
\examples{
dat <- cbind(rnorm(100),rep(1,100))
colnames(dat) <- c("y","x")
dat[42,1] <- NA
dat[57,2] <- NA
dat2 <- mxGREMLDataHandler(data=dat, yvars="y", Xvars=list("x"),
  addOnes = FALSE)
str(dat2)
}
