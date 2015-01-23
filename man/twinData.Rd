%
%   Copyright 2007-2015 The OpenMx Project
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

\name{twinData}
\alias{twinData}
\title{Australian twin sample biometric data.}
\description{
Australian twin data with 3808 observations on the 12 variables including body mass index (BMI) assessed in both MZ and DZ twins.

Questionnaires were mailed to 5967 pairs age 18 years and over. These data consist of completed questionnaires returned by both members of 3808 (64 percent) pairs. There are two cohort blocks in the data: a younger group (zyg 1:5), and an older group (zyg 6:10)
 
It is a wide dataset, with two individuals per line. Data include zygosity (zyg), along with heights in metres, weights in kg, and the derived variables BMI in kg/m^2 (stored as \dQuote{htwt1} and \dQuote{htwt2}), as well as the log of this variable, stored here as bm1 and bm2. The logged values are more closely normally distributed.

fam is a family identifier. Age is entered only once, as the both twins in each pair share a common age.

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

\references{
Martin, N. G. & Jardine, R. (1986). Eysenck's contribution to behavior genetics. In S. Modgil & C. Modgil (Eds.), \emph{Hans Eysenck: Consensus and Controversy.} Falmer Press: Lewes, Sussex.

Martin, N. G., Eaves, L. J., Heath, A. C., Jardine, R., Feindgold, L. M., & Eysenck, H. J. (1986). Transmission of social attitudes. \emph{Proceedings of the National Academy of Science}, \bold{83}, 4364-4368.
}

\usage{data(twinData)}
\format{
 A data frame with 3808 observations on the following 12 variables.
 \describe{
 \item{\code{fam}}{a numeric vector} of family IDs
 \item{\code{age}}{a numeric vector} of ages (years)
 \item{\code{zyg}}{a numeric vector} of zygosity (see below for important details)
 \item{\code{part}}{a numeric vector}
 \item{\code{wt1}}{a numeric vector} of weights in kg (twin 1)
 \item{\code{wt2}}{a numeric vector} of weights in kg (twin 2)
 \item{\code{ht1}}{a numeric vector} of heights in kg (twin 1)
 \item{\code{ht2}}{a numeric vector} of heights in kg (twin 2)
 \item{\code{htwt1}}{a numeric vector} of kg/m^2 twin 1
 \item{\code{htwt2}}{a numeric vector} of kg/m^2 twin 2
 \item{\code{bmi1}}{a numeric vector} of log BMI for twin 1
 \item{\code{bmi2}}{a numeric vector} of log BMI for twin 2
 }
}
\details{
Zygosity is coded as follows
1 == MZFF (i.e MZ females)
2 == MZMM (i.e MZ males)
3 == DZFF
4 == DZMM
5 == DZOS opposite sex pairs

Note: Zygosity 6:10 is the same, for an older cohort in the sample. So:
6 == MZFF (i.e MZ females)
7 == MZMM (i.e MZ males)
8 == DZFF
9 == DZMM
10 == DZOS opposite sex pairs

}
\examples{
data(twinData)
str(twinData)
plot(wt1 ~ wt2, data = twinData)
mzData <- as.matrix(subset(myTwinData, zyg == 1, c(bmi1, bmi2)))
dzData <- as.matrix(subset(myTwinData, zyg == 3, c(bmi1, bmi2)))
}
\keyword{datasets}
