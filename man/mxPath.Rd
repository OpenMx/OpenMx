\name{mxPath}
\alias{mxPath}
\alias{MxPath}

\title{Not Supported:Function To Create List of Paths}

\description{
   This function creates a list of paths.
}

\usage{
mxPath(from, to = NULL, all = FALSE, free = NULL, 
	arrows = NULL, startVal = NULL, 
	endVal = NULL, algebra = NULL,
	name = NULL, label = NULL,
	boundMax = NULL, boundMin = NULL,
	ciUpper = NULL, ciLower = NULL)
}

\arguments{
   \item{from}{character vector. these are the sources of the new paths.}
   \item{to}{character vector. these are the sinks of the new paths.}
   \item{all}{boolean. If TRUE, then connect all sources to all sinks.}
   \item{free}{boolean. Indicates whether paths are free or fixed.}
   \item{arrows}{numeric. Must be either 1 for single-headed or 2 for double-headed arrows.}
   \item{startVal}{numeric. The starting values of the parameters.}
   \item{endVal}{numeric. TODO what does this do?}
   \item{algebra}{character vector. The names of algebras associated with the paths.}
   \item{name}{character vector. The names of the paths.}
   \item{label}{character vector. Descriptive labels for the paths.}
   \item{boundMax}{numeric. TODO what does this do?}
   \item{boundMin}{numeric. TODO what does this do?}
   \item{ciUpper}{numeric. TODO what does this do?}
   \item{ciLower}{numeric. TODO what does this do?}
}

\details{
Not supported. mxPath is currently functional, but does not currently work with the mxRAMObjective function. Models are estimated, but mxPath objects are not turned into free parameters, and the resulting output includes only numeric(0) for parameter estimates.
}

\value{
	Returns a list of paths.
}

\examples{

myManifest <- sprintf("\%02d", c(1:100))
myLatent <- c("G1", "G2", "G3", "G4", "G5")
model <- mxModel()

for (i in 1:5) {
    j <- i*20
	singles <- mxPath(from=myLatent[i], to=myManifest[(j-19):j], 
    	arrows=1, free=c(FALSE,rep(TRUE, 19)), 
    	startVal=c(1,rep(0.75,19)))
	model <- mxModel(model, singles)
}

doubles <- mxPath(from=myLatent, all=TRUE, arrows=2,
    free=TRUE, startVal=1)

model <- mxModel(model, doubles)
}
