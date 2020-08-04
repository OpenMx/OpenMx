\name{imxGetNumThreads}
\alias{imxGetNumThreads}
\title{imxGetNumThreads}
\usage{
imxGetNumThreads()
}
\description{
This is an internal function exported for those people who know
what they are doing.

This function hard codes responses to a set of environments, like detecting snowfall,
or running on a cluster where "OMP_NUM_THREADS" is set or otherwise returning 1 or 2 cores
to avoid consuming all the resources on CRAN's test machines during release cycles.

This makes it \emph{not} suitable for getting the number of available threads.

To get the number of cores available locally you want \code{\link{omxDetectCores}}
or perhaps the \code{detectCores} function in the parallel package.

}
