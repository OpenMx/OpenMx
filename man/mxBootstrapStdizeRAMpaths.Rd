\name{mxBootstrapStdizeRAMpaths}
\alias{mxBootstrapStdizeRAMpaths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap distribution of standardized RAM path coefficients
}
\description{
Uses the bootstrap distribution of a RAM model's free parameters to produce a bootstrap distribution of standardized path coefficients.  Model must have already been run through \code{\link{mxBootstrap}()}.
}
\usage{
mxBootstrapStdizeRAMpaths(model, bq=c(.25,.75), 
	method=c('bcbci','quantile'), returnRaw=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{An MxModel that uses \link[=mxExpectationRAM]{RAM expectation} and has already been run through \code{\link{mxBootstrap}()}.}
  \item{bq}{Quantiles corresponding to the lower and upper limits, respectively, of the bootstrap confidence interval.}
  \item{method}{Character string.  One of 'bcbci' or 'quantile'.}
  \item{returnRaw}{Logical; should the function return the raw bootstrapping results?  Defaults to \code{FALSE}, in which case a dataframe summarizing the results is returned.}
}
\details{
In effect, what \code{mxBootstrapStdizeRAMpaths()} does is use the point estimates from each bootstrap replication to produce \code{\link{mxStandardizeRAMpaths}()} output for each replication.  The output of \code{\link{mxStandardizeRAMpaths}()} has one entry for each nonzero path coefficient.  Therefore, it is possible (though unlikely) that the number of nonzero paths, or which elements of the \code{A} and \code{S} RAM matrices are nonzero, may vary among bootstrap replications.  Such an occurrence defies simple summary of the standardized paths' bootstrapping results.  In this case, a raw list of bootstrapping results is returned, with a warning, if \code{returnRaw=TRUE}, or an error is throw if \code{returnRaw=FALSE}.

\code{mxBootstrapStdizeRAMpaths()} ignores any submodels of \code{model}.  That is, it must be able to find, in the MxModel it is provided, a \link[=mxExpectationRAM]{RAM expectation} and an \link[=mxComputeBootstrap]{MxComputeBootstrap} object.    It can be run on submodels directly,
}
\value{
Under extraordinary circumstances described above, \code{mxBootstrapStdizeRAMpaths()} may return a list object.  In ordinary circumstances, if \code{returnRaw=FALSE} (default), it returns a dataframe containing, \emph{inter alia}, the standardized path coefficients as estimated from the real data, their bootstrap SEs, and the lower and upper limits of a bootstrap confidence interval.  If instead \code{returnRaw=TRUE}, \code{mxBootstrapStdizeRAMpaths()} returns a matrix containing the raw bootstrap results; this matrix has one column per nonzero path coefficient, and one row for each succesfully converged bootstrap replication.
}
\seealso{
\code{\link{mxBootstrap}()}, \code{\link{mxStandardizeRAMpaths}()}, \code{\link{mxBootstrapEval}}, \code{\link{mxSummary}}
}
\examples{
require(OpenMx)
data(myFADataRaw)
myFADataRaw <- myFADataRaw[,c("x1","x2","x3","x4","x5","x6")]
dataRaw      <- mxData( observed=myFADataRaw, type="raw" )
resVars      <- mxPath( from=c("x1","x2","x3","x4","x5","x6"), arrows=2,
                        free=TRUE, values=c(1,1,1,1,1,1),
                        labels=c("e1","e2","e3","e4","e5","e6") ) 
latVar       <- mxPath( from="F1", arrows=2,
                        free=TRUE, values=1, labels ="varF1" )
facLoads     <- mxPath( from="F1", to=c("x1","x2","x3","x4","x5","x6"), arrows=1,
                        free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE), values=c(1,1,1,1,1,1),
                        labels =c("l1","l2","l3","l4","l5","l6") )
means        <- mxPath( from="one", to=c("x1","x2","x3","x4","x5","x6","F1"), arrows=1,
                        free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE), 
                        values=c(1,1,1,1,1,1,0),
                        labels =c("meanx1","meanx2","meanx3",
                                  "meanx4","meanx5","meanx6",NA) ) 
oneFactorModel <- mxModel("Common Factor Model Path Specification", type="RAM",
                        manifestVars=c("x1","x2","x3","x4","x5","x6"), latentVars="F1",
                        dataRaw, resVars, latVar, facLoads, means)
oneFactorFit <- mxRun(oneFactorModel)    
set.seed(170505)

\dontrun{oneFactorBoot <- mxBootstrap(oneFactorFit)
mxBootstrapStdizeRAMpaths(oneFactorBoot)}
}

