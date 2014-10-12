\name{HS.ability.data}
\alias{HS.ability.data}
\alias{HS.fake.data}
\docType{data}

\title{
Holzinger and Swineford (1939) Ability data in 301 children from two schools
}
\description{
This classic data set contains of intelligence-test scores from 301 children on 26 distinct tests.
The data are also available in the MBESS package.

The tests cover mental speed, memory, mathematical-ability, spatial, and verbal ability as listed below.
}

\usage{data("HS.ability.data")}

\format{
  A data frame with 301 observations on the following 2 variables.
  \describe{
		\item{\code{id}}{ student ID number (int)}
		\item{\code{Gender}}{ Sex (Factor w/ 2 levels \dQuote{Female},\dQuote{Male}}
		\item{\code{grade}}{ Grade in school (integer 7 or 8)}
		\item{\code{agey}}{  Age in years (integer)}
		\item{\code{agem}}{  Age in months (integer)}
		\item{\code{school}}{ School attended (Factor w/2 levels \dQuote{Grant-White} and \dQuote{Pasteur})}

		\item{\code{addition}}{ A speed test (numeric)}
		\item{\code{code}}{ A speed test (numeric)}
		\item{\code{counting}}{ A speed test (numeric)}
		\item{\code{straight}}{ A speed test (numeric)}

		\item{\code{wordr}}{ A memory subtest}
		\item{\code{numberr}}{ A memory subtest}
		\item{\code{figurer}}{ A memory subtest}
		\item{\code{object}}{ A memory subtest}
		\item{\code{numberf}}{ A memory subtest}
		\item{\code{figurew}}{ A memory subtest}

		\item{\code{deduct}}{A mathematical subtest}
		\item{\code{numeric}}{A mathematical subtest}
		\item{\code{problemr}}{A mathematical subtest}
		\item{\code{series}}{A mathematical subtest}
		\item{\code{arithmet}}{A mathematical subtest}

		\item{\code{visual}}{ A spatial subtest}
		\item{\code{cubes}}{ A spatial subtest}
		\item{\code{paper}}{ A spatial subtest}
		\item{\code{flags}}{ A spatial subtest}
		\item{\code{paperrev}}{ A spatial subtest}
		\item{\code{flagssub}}{ A spatial subtest}

		\item{\code{general}}{ A verbal subtest}
		\item{\code{paragrap}}{ A verbal subtest}
		\item{\code{sentence}}{ A verbal subtest}
		\item{\code{wordc}}{ A verbal subtest}
		\item{\code{wordm}}{ A verbal subtest}
  }
}

\details{
The data are from children who differ in grade (seventh- and eighth-grade) and are nested in one of two schools (Pasteur and Grant-White). You will see it in use elsewhere, both in R (lavaan, MBESS), and in Joreskog (1969) reporting a cfa on the Grant-White school subject subset).

The last two tests are substitute versions for other tests. \code{paperrev} (a paper form board test) can substitute for \code{paper} and \code{flagssub} for the lozenges test \code{flags}.
}

\source{
Holzinger, K., and Swineford, F. (1939).
}

\references{
Holzinger, K., and Swineford, F. (1939). A study in factor analysis: The stability of a bifactor solution. \emph{Supplementary Educational Monograph}, no. \bold{48}. Chicago: University of Chicago Press.
Joreskog, K. G. (1969). A general approach to confirmatory maximum likelihood factor analysis. \emph{Psychometrika}, \bold{34}, 183-202.
}

\examples{
data(HS.ability.data)
str(HS.ability.data)
levels(HS.ability.data$school)
plot(flags ~ flagssub, data = HS.ability.data)
}

\keyword{datasets}
