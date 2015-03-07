\name{myTwinData}
\alias{myTwinData}
\docType{data}
\title{
Twin data on weight and height
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("myTwinData")}
\format{
  A data frame with 3808 observations on the following variables.
  \describe{
    \item{\code{fam}}{Family ID variable}
    \item{\code{age}}{Age of the twin pair.  Range: 17 to 88.}
    \item{\code{zyg}}{Integer codes for zygosity and gender combinations}
    \item{\code{part}}{}
    \item{\code{wt1}}{Weight in kilograms for twin 1}
    \item{\code{wt2}}{Weight in kilograms for twin 2}
    \item{\code{ht1}}{Height in meters for twin 1}
    \item{\code{ht2}}{Height in meters for twin 2}
    \item{\code{htwt1}}{Product of ht and wt for twin 1}
    \item{\code{htwt2}}{Product of ht and wt for twin 2}
    \item{\code{bmi1}}{Body Mass Index for twin 1}
    \item{\code{bmi2}}{Body Mass Index for twin 2}
  }
}

\details{
Height and weight are highly correlated, and each individually highly heritable.  These data present and opportunity for multivariate behavior genetics modeling.
}

\source{
Timothy Bates
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
data(myTwinData)

plot( ht1 ~ wt1, myTwinData)

}
\keyword{datasets}
