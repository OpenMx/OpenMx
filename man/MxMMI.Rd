\name{mxModelavgEstimates}
\alias{mxModelavgEstimates}
\alias{omxAkaikeWeights}

\title{
Information-Theoretic Model-Averaging and Multimodel Inference
}
\description{
\code{omxAkaikeWeights()} orders a list of \link[=MxModel]{MxModels} from best to worst AIC, reports their Akaike weights, and indicates which are in the confidence set for best-approximating model.  \code{mxModelavgEstimates()} calls \code{omxAkaikeWeights()} and includes its output, and also reports model-average point estimates and their standard errors.
}
\usage{
mxModelavgEstimates(reference=character(0), models=list(),
include=c("onlyFree","all"), refAsBlock=FALSE, covariances=list(), 
type=c("AIC","AICc"), confidLvl=0.95)

omxAkaikeWeights(models=list(), type=c("AIC","AICc"), confidLvl=0.95)
}

\arguments{
  \item{reference}{Vector of character strings referring to \link[=omxGetParameters]{parameters}, \link[=MxMatrix]{MxMatrices}, or \link[=MxAlgebra]{MxAlgebras} for which model-average estimates are to be computed. Defaults to \code{NULL}.  If a zero-length value is provided, only the output of \code{omxAkaikeWeights()} is returned, with a warning.}
  \item{models}{List of at least two \link[=MxModel]{MxModel} objects, each of which must be uniquely identified by the value of its \code{name} slot.  Defaults to an empty list.}
  \item{include}{Character string, either \code{"onlyFree"} (default) or \code{"all"}.  When calculating model-average estimates for a given reference quantity, should all the MxModels in \code{models} be included in the calculations, or only those in which the quantity is freely estimated?  See below, under "Details," for additional information.}
  \item{refAsBlock}{Logical; defaults to \code{FALSE}.  Function behavior when \code{TRUE} is Not Yet Implemented.}
  \item{covariances}{Optional list of repeated-sampling covariance matrices of free parameter estimates (possibly from bootstrapping or the sandwich estimator); defaults to an empty list.  A non-empty list must either be of the same length as \code{models}, or have named elements corresponding to MxModel names.  See below, under "Details," for additional information.}
  \item{type}{Character string specifying which information criterion to use: either \code{"AIC"} for the ordinary AIC (default), or \code{"AICc"} for Hurvich & Tsai's (1984) sample-size corrected AIC.}
  \item{confidLvl}{Numeric proportion specifying the desired coverage probability of the confidence set for best-approximating model among \code{models} (Burnham & Anderson, 2002).  Defaults to 0.95.}
}
\details{
If statistical inferences (hypothesis tests and confidence intervals) are the motivation for calculating model-average point estimates and their standard errors, then \code{include="onlyFree"} (the default) is recommended.  Note that, if models in which a quantity is held fixed are included in calculating the quantity's model-average estimate, then that estimate cannot even asymptotically be normally distributed.

If argument \code{covariances} is non-empty, then either it must be of the same length as argument \code{models}, or all of its elements must be named after an MxModel in \code{models} (an MxModel's name is the character string in its \code{name} slot).  If \code{covariances} is of the same length as \code{models} but lacks element names, \code{mxModelavgEstimates()} will assume that they are ordered so that the first element of \code{covariances} is to be used with the first MxModel, the second element is to be used with the second MxModel, and so on.  Otherwise, \code{mxModelavgEstimates()} assigns the elements of \code{covariances} to the MxModels by matching element names to model names.  If \code{covariances} doesn't provide a covariance matrix for a given MxModel--perhaps because it is empty, or only provides matrices for a nonempty proper subset of the MxModels--\code{mxModelavgEstimates()} will fall back onto its default behavior of calculating a covariance matrix from the Hessian matrix in the MxModel's output slot.

The matrices in \code{covariances} must have complete row and column names, equal to the free parameter labels of the corresponding MxModel.  These names indicate to which free parameter a given row or column corresponds.
}
\value{
\code{omxAkaikeWeights()} returns a dataframe, with one row for each element of \code{models}.  The rows are sorted by their MxModel's AIC (or AICc), from best to worst.  The dataframe has five columns:
\enumerate{
	\item \code{"model"}: Character string. The name of the MxModel.
	\item \code{"AIC"} or \code{"AICc"}: Numeric.  The MxModel's AIC or AICc.
	\item \code{"delta"}: Numeric.  The MxModel's AIC(c) minus the best (smallest) AIC(c) in the set of MxModels.
	\item \code{"AkaikeWeight"}: Numeric.  The MxModel's Akaike weight.  This column will sum to unity.
	\item \code{"inConfidenceSet"}: Character.  Will contain an asterisk if the MxModel is in the confidence set for best-approximating model.
}
The dataframe also has an attribute, \code{"unsortedModelNames"}, which contains the names of the MxModels in the same order as they appear in \code{models} (i.e., without sorting them by their AIC). 

If a zero-length value is provided for argument \code{reference}, then \code{mxModelavgEstimates()} returns only the output of \code{omxAkaikeWeights()}, with a warning.  Otherwise, \code{mxModelavgEstimates()} returns a list with four elements:
\enumerate{
	\item A numeric matrix with one row for each distinct quantity specified by \code{reference}, and two columns.  Its rows are named for the corresponding reference quantities.  Its first column, \code{"Estimate"}, contains the model-average point estimates.  Its second column, \code{"SE"}, contains the "model-unconditional" standard errors of the model-average point estimates.
	\item A numeric matrix with one row for each distinct quantity specified by \code{reference} (indicated by row name), and one column for each MxModel (indicated by column name).  Each element is an estimate of the given quantity, from the given MxModel.  Quantities that cannot be evaluated for a given MxModel are reported as \code{NA}.
	\item A numeric matrix just like the one in list element 2, except that its elements are the esimated sampling variances of the corresponding model-conditional point estimates in list element 2.  Variances for fixed quantities are reported as 0 if \code{include="all"}, and as \code{NA} if \code{include="onlyFree"}.
	\item The output from \code{omxAkaikeWeights()} .
}

}
\references{
TODO
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
\dontrun{
mxModelavgEstimates(
	reference=c("A","Udiag"), include="all",
	models=list(factorFit1,factorFit2))
}
}

