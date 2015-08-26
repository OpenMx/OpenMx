\name{LongitudinalOverdispersedCounts}
\alias{LongitudinalOverdispersedCounts}
\docType{data}
\title{
Longitudinal, Overdispersed Count Data
}
\description{
Four-timepoint longitudinal data generated from an arbitrary Monte Carlo simulation, for 1000 simulees.  The response variable is a discrete count variable.  There are three time-invariant covariates.  The data are available in both "wide" and "long" format.
}
\usage{data("LongitudinalOverdispersedCounts")}
\format{
  The "long" format dataframe, \code{longData}, has 4000 rows and the following variables (columns):
  \enumerate{
  	\item{\code{id}: Factor; simulee ID code.}
  	\item{\code{tiem}: Numeric; represents the time metric, wave of assessment.}
  	\item{\code{x1}: Numeric; time-invariant covariate.}
  	\item{\code{x2}: Numeric; time-invariant covariate.}
  	\item{\code{x3}: Numeric; time-invariant covariate.}
  	\item{\code{y}: Numeric; the response ("dependent") variable.}
  }
  The "wide" format dataset, \code{wideData}, is a numeric 1000x12 matrix containing the following variables (columns):
  \enumerate{
  	\item{\code{id}: Simulee ID code.}
  	\item{\code{x1}: Time-invariant covariate.}
  	\item{\code{x3}: Time-invariant covariate.}
  	\item{\code{x3}: Time-invariant covariate.}
  	\item{\code{y0}: Response at initial wave of assessment.}
  	\item{\code{y1}: Response at first follow-up.}
  	\item{\code{y2}: Response at second follow-up.}
  	\item{\code{y3}: Response at third follow-up.}
  	\item{\code{t0}: Time variable at initial wave of assessment (in this case, 0).}
  	\item{\code{t1}: Time variable at first follow-up (in this case, 1).}
  	\item{\code{t2}: Time variable at second follow-up (in this case, 2).}
  	\item{\code{t3}: Time variable at third follow-up (in this case, 3).}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(LongitudinalOverdispersedCounts)
head(wideData)
str(longData)
#Let's try ordinary least-squares (OLS) regression:
olsmod <- lm(y~tiem+x1+x2+x3, data=longData)
#We will see in the diagnostic plots that the residuals are poorly approximated by normality, and are 
#heteroskedastic.  We also know that the residuals are not independent of one another, because we have 
#repeated-measures data:
plot(olsmod)
#In the summary, it looks like all of the regression coefficients are significantly different from zero, but we
#know that because the assumptions of OLS regression are violated that we should not trust its results:
summary(olsmod)

#Let's try a generalized linear model (GLM).  We'll use the quasi-Poisson quasilikelihood function to see
#how well the y variable is approximated by a Poisson distribution (conditional on time and covariates):
glm.mod <- glm(y~tiem+x1+x2+x3, data=longData, family="quasipoisson")
#The estimate of the dispersion parameter should be about 1.0 if the data are conditionally Poisson.  
#We can see that it is actually greater than 2, indicating overdispersion:
summary(glm.mod)
}
\keyword{datasets}
