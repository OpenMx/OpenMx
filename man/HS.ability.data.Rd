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
  A data frame comprising 301 observations on 22 variables:
  \describe{
		\item{\code{id}}{ student ID number (int)}
		\item{\code{Gender}}{ Sex (Factor w/ 2 levels \dQuote{Female},\dQuote{Male}}
		\item{\code{grade}}{ Grade in school (integer 7 or 8)}
		\item{\code{agey}}{  Age in years (integer)}
		\item{\code{agem}}{  Age in months (integer)}
		\item{\code{school}}{ School attended (Factor w/2 levels \dQuote{Grant-White} and \dQuote{Pasteur})}

		\item{\code{addition}}{ A speed test of addition (numeric)}
		\item{\code{code}}{ A speed test (numeric)}
		\item{\code{counting}}{ A speed test of counting groups of dots (numeric)}
		\item{\code{straight}}{ A speed test discriminating straight and curved capitals (numeric)}

		\item{\code{wordr}}{ A memory subtest of word recognition}
		\item{\code{numberr}}{ A memory subtest of number recognition}
		\item{\code{figurer}}{ A memory subtest of figure recognition}
		\item{\code{object}}{ A memory subtest: object-number test}
		\item{\code{numberf}}{ A memory subtest: number-figure test}
		\item{\code{figurew}}{ A memory subtest: figure-word test}

		\item{\code{deduct}}{A mathematical subtest of deduction}
		\item{\code{numeric}}{A mathematical subtest of numerical puzzles}
		\item{\code{problemr}}{A mathematical subtest of problem reasoning}
		\item{\code{series}}{A mathematical subtest of series completion}
		\item{\code{arithmet}}{A mathematical subtest: Woody-McCall mixed fundamentals, form I test, test 24}

		\item{\code{visual}}{ A spatial subtest of visual perception}
		\item{\code{cubes}}{ A spatial subtest}
		\item{\code{paper}}{ A spatial subtest paper form board}
		\item{\code{flags}}{ A spatial subtest (also known as lozenges)}
		\item{\code{paperrev}}{ A spatial subtest additional paper form board test (can substitute for paper)}
		\item{\code{flagssub}}{ A spatial subtest additional lozenges test (can substitute for flags)}

		\item{\code{general}}{ A verbal subtest of general information}
		\item{\code{paragrap}}{ A verbal subtest of paragraph comprehension}
		\item{\code{sentence}}{ A verbal subtest of sentence completion}
		\item{\code{wordc}}{ A verbal subtest of word classification}
		\item{\code{wordm}}{ A verbal subtest of word meaning}
  }
}

\details{
The data are from children who differ in grade (seventh- and eighth-grade) and are nested in one of two schools (Pasteur and Grant-White). You will see it in use elsewhere, both in R (lavaan, MBESS), and in Joreskog (1969) reporting a CFA on the Grant-White school subject subset).

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
