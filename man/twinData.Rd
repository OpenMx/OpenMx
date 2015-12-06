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
Australian twin data with 3,808 observations on the 12 variables including body mass index (BMI) assessed in both MZ and DZ twins.

Questionnaires were mailed to 5,967 pairs age 18 years and over. These data consist of completed questionnaires returned by both members of 3,808 (64 percent) pairs. There are two cohort blocks in the data: a younger group (zyg 1:5), and an older group (zyg 6:10)
 
It is a wide dataset, with two individuals per line. Families are identified by the variable \dQuote{fam}.

Data include zygosity (zyg), along with heights in metres, weights in kg, and the derived variables BMI in kg/m^2 (stored as \dQuote{htwt1} and \dQuote{htwt2}), as well as the log of this variable, stored here as bm1 and bm2. The logged values are more closely normally distributed.

For convenience, zyg is broken out into separate \dQuote{zygosity} and \dQuote{cohort} factors. \dQuote{zygosity} is coded as a 5-level factor.
}

\references{
Martin, N. G. & Jardine, R. (1986). Eysenck's contribution to behavior genetics. In S. Modgil & C. Modgil (Eds.), \emph{Hans Eysenck: Consensus and Controversy.} Falmer Press: Lewes, Sussex.

Martin, N. G., Eaves, L. J., Heath, A. C., Jardine, R., Feindgold, L. M., & Eysenck, H. J. (1986). Transmission of social attitudes. \emph{Proceedings of the National Academy of Science}, \bold{83}, 4364-4368.
}

\usage{data(twinData)}
\format{
 A data frame with 3808 observations on the following 12 variables.
 \describe{
 	\item{\code{fam}}{The family ID}
    \item{\code{age}}{Age in years (of both twins)}
    \item{\code{zyg}}{Code for zygosity and cohort (see details)}
    \item{\code{part}}{A numeric vector}
    \item{\code{wt1}}{Weight of twin 1 (kg)}
    \item{\code{wt2}}{Weight of twin 2 (kg)}
    \item{\code{ht1}}{Height of twin 1 (m)}
    \item{\code{ht2}}{Height of twin 2 (m)}
    \item{\code{htwt1}}{Raw BMI of twin 1 (kg/m^2)}
    \item{\code{htwt2}}{Raw BMI of twin 2 (kg/m^2)}
    \item{\code{bmi1}}{log(BMI) of twin 1}
    \item{\code{bmi2}}{log(BMI) of twin 2}
    \item{\code{cohort}}{Either \dQuote{younger} or \dQuote{older}}
    \item{\code{zygosity}}{Zygosity factor with levels: MZFF, MZMM, DZFF, DZMM, DZOS}
    \item{\code{age1}}{Age of Twin 1}
    \item{\code{age2}}{Age of Twin 2}
 }
}

\details{
\dQuote{zyg} codes twin-zygosity as follows:
1 == MZFF (i.e MZ females)
2 == MZMM (i.e MZ males)
3 == DZFF
4 == DZMM
5 == DZOS opposite sex pairs

Note: zyg 6:10 are for an older cohort in the sample. So:
6 == MZFF (i.e MZ females)
7 == MZMM (i.e MZ males)
8 == DZFF
9 == DZMM
10 == DZOS opposite sex pairs

The \dQuote{zygosity} and \dQuote{cohort} variables take care of this for you (conventions differ).
}
\examples{
data(twinData)
str(twinData)
plot(wt1 ~ wt2, data = twinData)
selVars = c("bmi1", "bmi2")
mzData <- subset(twinData, zyg == 1, selVars)
dzData <- subset(twinData, zyg == 3, selVars)

# equivalently
mzData <- subset(twinData, zygosity == "MZFF", selVars)

# Disregard sex, pick older cohort
mz <- subset(twinData, zygosity \%in\% c("MZFF","MZMM") & cohort == "older", selVars)

}
\keyword{datasets}
