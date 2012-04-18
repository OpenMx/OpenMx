\name{twinData}
\alias{twinData}
\title{Australian twin sample biometric data.}
\description{
Australian data on body mass index (BMI) assessed in both MZ and DZ twins, and saved in the text file twinData.txt. 

It is a wide dataset, with two individuals per line. Data include zygosity (zyg), along with heights in metres, weights in kg, and the derived variables BMI in kg/m^2 (stored as "htwt1" and "htwt2"), as well as the log of this variable, stored here as bm1 and bm2. The logged values are more closely normally distributed.

Age of course occurs only once, as the two twins of each pair share a common age. fam is a family identifier.

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
Zygosity is coded as follows
1 == MZ females
2 == MZ males
3 == DZ females
4 == DZ males
5 == DZ opposite sex pairs
}
\examples{
data(twinData)
str(twinData)
plot(wt1~wt2, data=twinData)
mzData <- as.matrix(subset(myTwinData, zyg==1, c(bmi1,bmi2)))
dzData <- as.matrix(subset(myTwinData, zyg==3, c(bmi1,bmi2)))
}
\keyword{datasets}