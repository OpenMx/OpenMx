\name{nuclear_twin_design_data}
\alias{nuclear_twin_design_data}
\docType{data}
\title{
Twin data from a nuclear family design
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("nuclear_twin_design_data")}
\format{
  A data frame with 1743 observations on the following variables.
  \describe{
    \item{\code{Twin1}}{}
    \item{\code{Twin2}}{}
    \item{\code{Father}}{}
    \item{\code{Mother}}{}
    \item{\code{zyg}}{Zygosity of the twin pair}
  }
}

\details{
This is a wide format data set.  A single variable has values for different member of the same nuclear family.
}

\source{
Likely simulated.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(nuclear_twin_design_data)

cor(nuclear_twin_design_data[,-5], use="pairwise.complete.obs")

}
\keyword{datasets}
