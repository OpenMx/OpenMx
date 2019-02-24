\name{twin_NA_dot}
\alias{twin_NA_dot}
\docType{data}
\title{
Twin biometric data (Practice cleaning: "." for missing data, wrong data types etc.)
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("twin_NA_dot")}
\format{
  A data frame with 3808 observations on the following variables.
  \describe{
    \item{\code{fam}}{Family ID variable}
    \item{\code{age}}{Age of the twin pair.  Range: 17 to 88, coded as factor}
    \item{\code{zyg}}{Integer codes for zygosity and gender combinations}
    \item{\code{part}}{Cohort}
    \item{\code{wt1}}{Weight in kilograms for twin 1 (this and following have "." embedded as NA...)}
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
Same as \link{myTwinData} but has . as the missing data value instead of NA. 
}

\source{
Timothy Bates
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\examples{
data(twin_NA_dot)
summary(twin_NA_dot)
# Note that all variables are treated as factors because of the missing data coding.


}
\keyword{datasets}
