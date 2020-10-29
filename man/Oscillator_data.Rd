\name{Oscillator}
\alias{Oscillator}
\docType{data}
\title{
Oscillator Data for Latent Differential Equations
}

\description{
Data set used in some of OpenMx's examples, for instance the LDE demo. The data were simulated by Steven M. Boker according to a noisy oscillator model.
}

\usage{data("Oscillator")}
\format{
  A data frame with 100 observations on the following 1 numeric variable.
  \describe{
    \item{\code{x}}{Noisy oscillator value}
  }
}

\details{
The data appear to be sinusoidal with exponential decay on the amplitude.  The rows of data are different times.  The column is the variable.}

\source{
Simulated.  Pulled from \url{https://openmx.ssri.psu.edu/thread/144}}

\references{

Boker, S., Neale, M., & Rausch, J. (2004). Latent differential equation modeling with multivariate multi-occasion indicators. In \emph{Recent developments on structural equation models} van Montfort, K., Oud, H, and Satorra, A. (Eds.). 151-174. Springer, Dordrecht.

}

\examples{
data(Oscillator)
plot(Oscillator$x, type='l')
}
\keyword{datasets}
