\name{mxRestore}
\alias{mxRestore}

\title{Restore From Checkpoint File}

\description{
   The function loads the last saved state from a checkpoint file.
}

\usage{
mxRestore(model, chkpt.directory = ".", chkpt.prefix = "")
}

\arguments{
    \item{model}{\link{MxModel} object to be loaded.}    
	\item{chkpt.directory}{character. Directory where the checkpoint file is located.}
	\item{chkpt.prefix}{character. Prefix of the checkpoint file.}
}

\details{
In general, the arguments \sQuote{chkpt.directory} and \sQuote{chkpt.prefix} should be identical to the \code{\link{mxOption}}: \sQuote{Checkpoint Directory} and \sQuote{Checkpoint Prefix} that were specificed on the model before execution.

Alternatively, the checkpoint file can be manually loaded as a data.frame in R.  Use \code{\link{read.table}} with the options \sQuote{header=TRUE}, \sQuote{stringsAsFactors=FALSE} and \sQuote{check.names=FALSE}.
}

\value{
   Returns an MxModel object with free parameters updated to the last saved values.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
#Create a model that includes data, 
#matrices A, S and F, and an objective function
\dontrun{
data <- mxData(mydata, type="cov", numObs = 100)
objective <- mxRAMObjective('A', 'S', 'F')
model <- mxModel("mymodel", A, S, F, data, objective)

#Use mxRun to optimize the free parameters in the matrices A and S
modelOut <- mxRun(model, checkpoint = TRUE)

#Use mxRestore to load the last saved state of the model
modelRestore <- mxRestore(model)
}
}
