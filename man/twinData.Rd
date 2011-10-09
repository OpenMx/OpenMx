\name{twinData}
\alias{twinData}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Australian twin sample biometric data.
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Australian data on body mass index (BMI) which are saved in a text file ‘twinData.txt’. It is a "wide" dataset, with two individuals per line. It also contains both MZ and DZ twins, along with heights, weights, and the calculated variable body mass index (BMI) for each subject. Age of course occurs only once, as the two twins of each pair share a common age. fam is a family identifier.

}
\usage{data(twinData)}
\format{
  A data frame with 3808 observations on the following 12 variables.
  \describe{
    \item{\code{fam}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{zyg}}{a numeric vector}
    \item{\code{part}}{a numeric vector}
    \item{\code{wt1}}{a numeric vector}
    \item{\code{wt2}}{a numeric vector}
    \item{\code{ht1}}{a numeric vector}
    \item{\code{ht2}}{a numeric vector}
    \item{\code{htwt1}}{a numeric vector}
    \item{\code{htwt2}}{a numeric vector}
    \item{\code{bmi1}}{a numeric vector}
    \item{\code{bmi2}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
Zygosity is coded as follows:
1 == MZ females
2 == MZ males
3 == DZ females
4 == DZ males
4 == DZ opposite sex pairs

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(twinData)
str(twinData) ; plot(wt1, wt2, data=twinData)
mzData <- as.matrix(subset(myTwinData, zyg==1, c(bmi1,bmi2)))
dzData <- as.matrix(subset(myTwinData, zyg==3, c(bmi1,bmi2)))
}
\keyword{datasets}
