\name{mxModelAverage}
\alias{mxModelAverage}
\alias{omxAkaikeWeights}
\alias{omxAICWeights}

\title{
Information-Theoretic Model-Averaging and Multimodel Inference
}
\description{
\code{omxAkaikeWeights()} orders a list of \link[=MxModel]{MxModels} (hereinafter, the "candidate set" of models) from best to worst AIC, reports their Akaike weights, and indicates which are in the confidence set for best-approximating model.  \code{mxModelAverage()} calls \code{omxAkaikeWeights()} and includes its output, and also reports model-average point estimates and (if requested) their standard errors.
}
\usage{
mxModelAverage(reference=character(0), models=list(),
include=c("onlyFree","all"), SE=NULL, refAsBlock=FALSE, covariances=list(), 
type=c("AIC","AICc"), conf.level=0.95)

omxAkaikeWeights(models=list(), type=c("AIC","AICc"), conf.level=0.95)
}

\arguments{
  \item{reference}{Vector of character strings referring to \link[=omxGetParameters]{parameters}, \link[=MxMatrix]{MxMatrices}, or \link[=MxAlgebra]{MxAlgebras} for which model-average estimates are to be computed. Defaults to \code{NULL}.  If a zero-length value is provided, only the output of \code{omxAkaikeWeights()} is returned, with a warning.}
  \item{models}{The candidate set of models: a list of at least two \link[=MxModel]{MxModel} objects, each of which must be uniquely identified by the value of its \code{name} slot.  Defaults to an empty list.}
  \item{include}{Character string, either \code{"onlyFree"} (default) or \code{"all"}.  When calculating model-average estimates for a given reference quantity, should all the MxModels in the candidate set be included in the calculations, or only those in which the quantity is freely estimated?  See below, under "Details," for additional information.}
  \item{SE}{Logical; should standard errors be reported for the model-average point estimates?  Defaults to \code{NULL}, in which case standard errors are reported if argument \code{include="onlyFree"}, and not reported otherwise.}
  \item{refAsBlock}{Logical. If \code{FALSE} (default), \code{mxModelAverage()} will include a matrix of model-conditional sampling variances for the reference quantities in its output, and model-average results may be based on different subsets of the candidate set if \code{include="onlyFree"}.  If \code{TRUE}, \code{mxModelAverage()} will instead include a joint sampling covariance matrix for all reference quantities, and will throw an error if \code{include="onlyFree"} and if it is not the case that all reference quantities are freely estimated in all models in the candidate set.}
  \item{covariances}{Optional list of repeated-sampling covariance matrices of free parameter estimates (possibly from bootstrapping or the sandwich estimator); defaults to an empty list.  A non-empty list must either be of the same length as \code{models}, or have named elements corresponding to names of MxModels in the candidate set.  See below, under "Details," for additional information.}
  \item{type}{Character string specifying which information criterion to use: either \code{"AIC"} for the ordinary AIC (default), or \code{"AICc"} for Hurvich & Tsai's (1989) sample-size corrected AIC.}
  \item{conf.level}{Numeric proportion specifying the desired coverage probability of the confidence set for best-approximating model among the candidate set (Burnham & Anderson, 2002).  Defaults to 0.95.}
}

\details{
If statistical inferences (hypothesis tests and confidence intervals) are the motivation for calculating model-average point estimates and their standard errors, then \code{include="onlyFree"} (the default) is recommended.  Note that, if models in which a quantity is held fixed are included in calculating the quantity's model-average estimate, then that estimate cannot even asymptotically be normally distributed (Bartels, 1997).

If argument \code{covariances} is non-empty, then either it must be of the same length as argument \code{models}, or all of its elements must be named after an MxModel in \code{models} (an MxModel's name is the character string in its \code{name} slot).  If \code{covariances} is of the same length as \code{models} but lacks element names, \code{mxModelAverage()} will assume that they are ordered so that the first element of \code{covariances} is to be used with the first MxModel, the second element is to be used with the second MxModel, and so on.  Otherwise, \code{mxModelAverage()} assigns the elements of \code{covariances} to the MxModels by matching element names to MxModel names.  If \code{covariances} doesn't provide a covariance matrix for a given MxModel--perhaps because it is empty, or only provides matrices for a nonempty proper subset of the candidate set--\code{mxModelAverage()} will fall back to its default behavior of calculating a covariance matrix from the Hessian matrix in the MxModel's output slot.  If a covariance matrix cannot be thus calculated and \code{SE=TRUE}, \code{SE} is coerced to \code{FALSE}, with a warning.

The matrices in \code{covariances} must have complete row and column names, equal to the free parameter labels of the corresponding MxModel.  These names indicate to which free parameter a given row or column corresponds.
}

\value{
\code{omxAkaikeWeights()} returns a dataframe, with one row for each element of \code{models}.  The rows are sorted by their MxModel's AIC (or AICc), from best to worst.  The dataframe has five columns:
\enumerate{
	\item \code{"model"}: Character string. The name of the MxModel.
	\item \code{"AIC"} or \code{"AICc"}: Numeric.  The MxModel's AIC or AICc.
	\item \code{"delta"}: Numeric.  The MxModel's AIC (or AICc) minus the best (smallest) AIC (or AICc) in the candidate set.
	\item \code{"AkaikeWeight"}: Numeric.  The MxModel's Akaike weight.  This column will sum to unity.
	\item \code{"inConfidenceSet"}: Character.  Will contain an asterisk if the MxModel is in the confidence set for best-approximating model.
}
The dataframe also has an attribute, \code{"unsortedModelNames"}, which contains the names of the MxModels in the same order as they appear in \code{models} (i.e., without sorting them by their AIC). 

If a zero-length value is provided for argument \code{reference}, then \code{mxModelAverage()} returns only the output of \code{omxAkaikeWeights()}, with a warning.  Otherwise, for the default values of its arguments, \code{mxModelAverage()} returns a list with four elements:
\enumerate{
	\item \code{"Model-Average Estimates"}: A numeric matrix with one row for each distinct quantity specified by \code{reference}, and as many as two columns.  Its rows are named for the corresponding reference quantities.  Its first column, \code{"Estimate"}, contains the model-average point estimates.  If standard errors are being calculated, then its second column, \code{"SE"}, contains the "model-unconditional" standard errors of the model-average point estimates.  Otherwise, there is no second column.
	\item \code{"Model-wise Estimates"}: A numeric matrix with one row for each distinct quantity specified by \code{reference} (indicated by row name), and one column for each MxModel (indicated by column name).  Each element is an estimate of the given reference quantity, from the given MxModel.  Quantities that cannot be evaluated for a given MxModel are reported as \code{NA}.
	\item \code{"Model-wise Sampling Variances"}: A numeric matrix just like the one in list element 2, except that its elements are the estimated sampling variances of the corresponding model-conditional point estimates in list element 2.  Variances for fixed quantities are reported as 0 if \code{include="all"}, and as \code{NA} if \code{include="onlyFree"}; however, if no covariance matrix is available for a model, all of that model's sampling variances will be reported as \code{NA}.
	\item \code{"Akaike-Weights Table"}: The output from \code{omxAkaikeWeights()}.
}
If \code{refAsBlock=TRUE}, list element 3 will instead contain be named \code{"Joint Covariance Matrix"}, and if \code{SE=TRUE}, it will contain the joint sampling covariance matrix for the model-average point estimates.
}

\note{
The "best-approximating model" is defined as the model that truly ("in the population," so to speak) has the smallest Kullback-Leibler divergence from full reality, among the models in the candidate set (Burnham & Anderson, 2002).

A model's Akaike weight is interpretable as the relative weight-of-evidence for that model being the best-approximating model, given the observed data and the candidate set.  It has a Bayesian interpretation as the posterior probability that the given model is the best-approximating model in the candidate set, assuming a "savvy" prior probability that depends upon sample size and the number of free parameters in the model (Burnham & Anderson, 2002).

The confidence set for best-approximating model serves to reflect sampling error in the AICs.  When fitting the candidate set to data over repeated sampling, the confidence set is expected to contain the best-approximating model with probability equal to its confidence level.

The sampling variances and covariances of the model-average point estimates are calculated from Equations (4) and (5) in Burnham & Anderson (2004).  The standard errors reported by \code{mxModelAverage()} are the square roots of those sampling variances.

For an example of model-averaging and multimodel inference applied to structural equation modeling using OpenMx v1.3 (i.e., well before the functions documented here were implemented), see Kirkpatrick, McGue, & Iacono (2015).
}

\references{
Bartels, L. M. (1997).  Specification uncertainty and model averaging.  \emph{American Journal of Political Science, 41}(2), 641-674.

Burnham, K. P., & Anderson, D. R.  (2002).  \emph{Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach (2nd ed.)}.  New York: Springer.

Burnham, K. P., & Anderson, D. R.  (2004).  Multimodel inference: Understanding AIC and BIC in model selection.  \emph{Sociological Methods & Research, 33}(2), 261-304.  doi:10.1177/0049124104268644

Hurvich, C. M., & Tsai, C-L.  (1989).  Regression and time series model selection in small samples.  \emph{Biometrika, 76}(2), 297-307.

Kirkpatrick, R. M., McGue, M., & Iacono, W. G.  (2015).  Replication of a gene-environment interaction via multimodel inference: Additive-genetic variance in adolescents' general cognitive ability increases with family-of-origin socioeconomic status.  \emph{Behavior Genetics, 45}, 200-214.
}

\seealso{
\code{\link{mxCompare}()}
}
\examples{
require(OpenMx)
data(demoOneFactor)
factorModel1 <- mxModel(
	"OneFactor1",
	mxMatrix(
		"Full", 5, 1, values=0.8, 
		labels=paste("a",1:5,sep=""),
		free=TRUE, name="A"),
	mxMatrix(
		"Full", 5, 1, values=1,
		labels=paste("u",1:5,sep=""),
		free=TRUE, name="Udiag"),
	mxMatrix(
		"Symm", 1, 1, values=1,
		free=FALSE, name="L"),
	mxAlgebra(vec2diag(Udiag),name="U"),
	mxAlgebra(A \%*\% L \%*\% t(A) + U, name="R"),
	mxExpectationNormal(
		covariance = "R",
		dimnames = names(demoOneFactor)),
	mxFitFunctionML(),
	mxData(cov(demoOneFactor), type="cov", numObs=500))
factorFit1 <- mxRun(factorModel1)
#Constrain unique variances equal:
factorModel2 <- omxSetParameters(
	model=factorModel1,labels=paste("u",1:5,sep=""),
	newlabels="u",name="OneFactor2")
factorFit2 <- mxRun(factorModel2)
omxAkaikeWeights(models=list(factorFit1,factorFit2))
\donttest{
mxModelAverage(
	reference=c("A","Udiag"), include="all",
	models=list(factorFit1,factorFit2))
}
}

