\name{Bollen}
\alias{Bollen}
\docType{data}
\title{
Bollen Data on Industrialization and Political Democracy
}

\description{
Data set used in some of OpenMx's examples, for instance WLS. The data were reported in Bollen (1989, p. 428, Table 9.4)
This set includes data from 75 developing countries each assessed on four measures of democracy measured twice (1960 and 1965), and three measures of industrialization measured once (1960).
}

\usage{data("Bollen")}
\format{
  A data frame with 75 observations on the following 11 numeric variables.
  \describe{
    \item{\code{y1}}{Freedom of the press, 1960 }
    \item{\code{y2}}{Freedom of political opposition, 1960}
    \item{\code{y3}}{Fairness of elections, 1960}
    \item{\code{y4}}{Effectiveness of elected legislature, 1960}
    \item{\code{y5}}{Freedom of the press, 1965}
    \item{\code{y6}}{Freedom of political opposition, 1965}
    \item{\code{y7}}{Fairness of elections, 1965}
    \item{\code{y8}}{Effectiveness of elected legislature, 1965}
    \item{\code{x1}}{GNP per capita, 1960}
    \item{\code{x2}}{Energy consumption per capita, 1960}
    \item{\code{x3}}{Percentage of labor force in industry, 1960}
  }
}

\details{
Variables y1-y4 and y5-y8 are typically used as indicators of the latent trait of \dQuote{political democracy} in 1960 and 1965 respectively. x1-x3 are used as indicators of industrialization (1960).
}

\source{
The \code{sem} package (in turn, via pers. comm Bollen to Fox)}

\references{

Bollen, K. A. (1979). Political democracy and the timing of development. \emph{American Sociological Review}, \bold{44}, 572-587.

Bollen, K. A. (1980). Issues in the comparative measurement of political democracy. \emph{American Sociological Review}, \bold{45}, 370-390.

Bollen, K. A. (1989). \emph{Structural equation models}. New York: Wiley-Interscience.

}

\examples{
data(Bollen)
str(Bollen)
plot(y1 ~ y2, data = Bollen)
}
\keyword{datasets}