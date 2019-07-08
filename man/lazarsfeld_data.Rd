\name{lazarsfeld}
\alias{lazarsfeld}
\docType{data}
\title{
Respondent-soldiers on four dichotomous items
}

\description{
Data set used in some of OpenMx's examples.
}

\usage{data("lazarsfeld")}
\format{
  A data frame with 1000 observations on four dichotomous items.
  \describe{
    \item{\code{armyrun}}{In general how do you feel the Army is run?}
    \item{\code{favatt}}{Do you think when you are discharged you will [have] a favorable attitude toward the Army?}
    \item{\code{squaredeal}}{In general do you feel you yourself have gotten a square deal from the Army?}
    \item{\code{welfare}}{Do you feel that the Army is trying its best to look out for the welfare of enlisted men?}
    \item{\code{frequency}}{Frequency of response pattern.}
  }
}

\details{
A straightforward descriptive analysis of these data shows that negative responses are more numerous except on item 1; and that there is a positive association between each pair of items. A soldier who responds positively to any one item is more likely to respond positively to a second item. Lazarsfeld's analysis is based on the assumption that each soldier can be thought of as belong to one of two latent classes. The probability of positive response to an item is different in one group than in the other. Most importantly, he is willing to assume that for an individual respondent the responses to items are statistically independent.
}

\source{
Lazarsfeld, Paul F. (1950b) "Some Latent Structures", Chapter 11 in Stouffer (1950).
}

\seealso{
http://www.people.vcu.edu/~nhenry/LSA50.htm
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\examples{
data(lazarsfeld)

}
\keyword{datasets}
