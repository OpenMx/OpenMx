%
%   Copyright 2007-2014 The OpenMx Project
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

\name{mxExpectationLISREL}
\alias{mxExpectationLISREL}
\alias{MxExpectationLISREL-class}
\alias{show,MxExpectationLISREL-method}
\alias{print,MxExpectationLISREL-method}

\title{Create MxExpectationLISREL Object}

\description{
   This function creates a new MxExpectationLISREL object.
}

\usage{
mxExpectationLISREL(LX=NA, LY=NA, BE=NA, GA=NA, PH=NA, PS=NA, TD=NA, TE=NA, TH=NA,
                    TX = NA, TY = NA, KA = NA, AL = NA,
                    dimnames = NA, thresholds = NA, threshnames = dimnames)
}


\arguments{
   \item{LX}{An optional character string indicating the name of the 'LX' matrix.}
   \item{LY}{An optional character string indicating the name of the 'LY' matrix.}
   \item{BE}{An optional character string indicating the name of the 'BE' matrix.}
   \item{GA}{An optional character string indicating the name of the 'GA' matrix.}
   \item{PH}{An optional character string indicating the name of the 'PH' matrix.}
   \item{PS}{An optional character string indicating the name of the 'PS' matrix.}
   \item{TD}{An optional character string indicating the name of the 'TD' matrix.}
   \item{TE}{An optional character string indicating the name of the 'TE' matrix.}
   \item{TH}{An optional character string indicating the name of the 'TH' matrix.}
   \item{TX}{An optional character string indicating the name of the 'TX' matrix.}
   \item{TY}{An optional character string indicating the name of the 'TY' matrix.}
   \item{KA}{An optional character string indicating the name of the 'KA' matrix.}
   \item{AL}{An optional character string indicating the name of the 'AL' matrix.}
   \item{dimnames}{An optional character vector that is currently ignored}%to be assigned to the column names of the 'F' and 'M' matrices.}   
   \item{thresholds}{An optional character string indicating the name of the thresholds matrix.}
   \item{threshnames}{An optional character vector to be assigned to the column names of the thresholds matrix.}
}

\details{
Expectation functions define the way that model expectations are calculated. The mxExpectationLISREL calculates the expected covariance and means of a given \link{MxData} object given a LISREL model. This model is defined by LInear Structual RELations (LISREL; J\ifelse{latex}{\out{\"o}}{\ifelse{html}{\out{&ouml;}}{o}}reskog & S\ifelse{latex}{\out{\"o}}{\ifelse{html}{\out{&ouml;}}{o}}rbom, 1982, 1996).  Arguments 'LX' through 'AL' must refer to \link{MxMatrix} objects with the associated properties of their respective matrices in the LISREL modeling approach.

The full LISREL specification has 13 matrices and is sometimes called the extended LISREL model.  It is defined by the following equations.

\deqn{\eta = \alpha + B \eta + \Gamma \xi + \zeta}
\deqn{y = \tau_y + \Lambda_y \eta + \epsilon}
\deqn{x = \tau_x + \Lambda_x \xi + \delta}

The table below is provided as a quick reference to the numerous matrices in LISREL models.  Note that NX is the number of manifest exogenous (independent) variables, the number of Xs.  NY is the number of manifest endogenous (dependent) variables, the number of Ys.  NK is the number of latent exogenous variables, the number of Ksis or Xis.  NE is the number of latent endogenous variables, the number of etas.

\tabular{clcccl}{
Matrix \tab Word \tab Abbreviation \tab Dimensions \tab Expression\tab Description \cr
\eqn{\Lambda_x} \tab Lambda x \tab LX \tab NX x NK \tab \tab Exogenous Factor Loading Matrix \cr
\eqn{\Lambda_y} \tab Lambda y \tab LY \tab NY x NE \tab \tab Endogenous Factor Loading Matrix\cr
\eqn{B} \tab Beta \tab BE \tab NE x NE \tab \tab Regressions of Latent Endogenous Variables Predicting Endogenous Variables\cr
\eqn{\Gamma} \tab Gamma \tab GA \tab NE x NK \tab \tab Regressions of Latent Exogenous Variables Predicting Endogenous Variables \cr
\eqn{\Phi} \tab Phi \tab PH \tab NK x NK \tab cov(\eqn{\xi}) \tab Covariance Matrix of Latent Exogenous Variables \cr
\eqn{\Psi} \tab Psi \tab PS \tab NE x NE \tab cov(\eqn{\zeta}) \tab Residual Covariance Matrix of Latent Endogenous Variables \cr
\eqn{\Theta_{\delta}} \tab Theta delta \tab TD \tab NX x NX \tab cov(\eqn{\delta}) \tab Residual Covariance Matrix of Manifest Exogenous Variables \cr
\eqn{\Theta_{\epsilon}} \tab Theta epsilon \tab TE \tab NY x NY \tab cov(\eqn{\epsilon}) \tab Residual Covariance Matrix of Manifest Endogenous Variables \cr
\eqn{\Theta_{\delta \epsilon}} \tab Theta delta epsilson \tab TH \tab NX x NY \tab cov(\eqn{\delta, \epsilon}) \tab Residual Covariance Matrix of Manifest Exogenous with Endogenous Variables \cr
\eqn{\tau_x} \tab tau x \tab TX \tab NX x 1 \tab \tab Residual Means of Manifest Exogenous Variables \cr
\eqn{\tau_y} \tab tau y \tab TY \tab NY x 1 \tab \tab Residual Means of Manifest Endogenous Variables \cr
\eqn{\kappa} \tab kappa \tab KA \tab NK x 1 \tab mean(\eqn{\xi}) \tab Means of Latent Exogenous Variables \cr
\eqn{\alpha} \tab alpha \tab AL \tab NE x 1 \tab \tab Residual Means of Latent Endogenous Variables \cr
}

From the extended LISREL model, several submodels can be defined.  Subtypes of the LISREL model are defined by setting some of the arguments of the LISREL expectation function to NA.  Note that because the default values of each LISREL matrix is NA, setting a matrix to NA can be accomplished by simply not giving it any other value.

The first submodel is the LISREL model without means.

\deqn{\eta = B \eta + \Gamma \xi + \zeta}
\deqn{y = \Lambda_y \eta + \epsilon}
\deqn{x = \Lambda_x \xi + \delta}

The LISREL model without means requires 9 matrices: LX, LY, BE, GA, PH, PS, TD, TE, and TH.  Hence this LISREL model has TX, TY, KA, and AL as NA.  This can be accomplished be leaving these matrices at their default values.

The TX, TY, KA, and AL matrices must be specified if either the mxData type is \dQuote{cov} or \dQuote{cor} and a means vector is provided, or if the mxData type is \dQuote{raw}.  Otherwise the TX, TY, KA, and AL matrices are ignored and the model without means is estimated.

A second submodel involves only endogenous variables.

\deqn{\eta = B \eta + \zeta}
\deqn{y = \Lambda_y \eta + \epsilon}

The endogenous-only LISREL model requires 4 matrices: LY, BE, PS, and TE.  The LX, GA, PH, TD, and TH must be NA in this case.  However, means can also be specified, allowing TY and AL if the data are raw or if observed means are provided.

Another submodel involves only exogenous variables.

\deqn{x = \Lambda_x \xi + \delta}

The exogenous-model model requires 3 matrices: LX, PH, and TD.  The LY, BE, GA, PS, TE, and TH matrices must be NA.  However, means can also be specified, allowing TX and KA if the data are raw or if observed means are provided.

The model that is run depends on the matrices that are not NA.  If all 9 matrices are not NA, then the full model is run.  If only the 4 endogenous matrices are not NA, then the endogenous-only model is run.  If only the 3 exogenous matrices are not NA, then the exogenous-only model is run.  If some endogenous and exogenous matrices are not NA, but not all of them, then appropriate errors are thrown.  Means are included in the model whenever their matrices are provided.

%The 'dimnames' arguments takes an optional character vector.  If this argument is not a single NA, then this vector be assigned to be the column names of the 'F' matrix and optionally to the 'M' matrix, if the 'M' matrix exists.

The \link{MxMatrix} objects included as arguments may be of any type, but should have the properties described above. The mxExpectationLISREL will not return an error for incorrect specification, but incorrect specification will likely lead to estimation problems or errors in the \link{mxRun} function.

Like the \link{mxExpectationRAM}, the mxExpectationLISREL evaluates with respect to an \link{MxData} object. The \link{MxData} object need not be referenced in the mxExpectationLISREL function, but must be included in the \link{MxModel} object. mxExpectationLISREL requires that the 'type' argument in the associated \link{MxData} object be equal to 'cov', 'cor', or 'raw'.

To evaluate, place mxExpectationLISREL objects, the \link{mxData} object for which the expected covariance approximates, referenced \link{MxAlgebra} and \link{MxMatrix} objects, and optional \link{MxBounds} and \link{MxConstraint} objects in an \link{MxModel} object. This model may then be evaluated using the \link{mxRun} function. The results of the optimization can be found in the 'output' slot of the resulting model, and may be obtained using the \link{mxEval} function.
}

\value{
    Returns a new MxExpectationLISREL object. One and only one MxExpectationLISREL object can be included with models using one and only one fit function object (e.g., MxFitFunctionML) and with referenced \link{MxAlgebra}, \link{MxData} and \link{MxMatrix} objects.
}

\references{
J\ifelse{latex}{\out{\"o}}{\ifelse{html}{\out{&ouml;}}{o}}reskog, K. G. & S\ifelse{latex}{\out{\"o}}{\ifelse{html}{\out{&ouml;}}{o}}rbom, D. (1996). LISREL 8: User's Reference Guide. Lincolnwood, IL: Scientific Software International.

J\ifelse{latex}{\out{\"o}}{\ifelse{html}{\out{&ouml;}}{o}}reskog, K. G. & S\ifelse{latex}{\out{\"o}}{\ifelse{html}{\out{&ouml;}}{o}}rbom, D. (1982). Recent developments in structural equation modeling. \emph{Journal of Marketing Research, 19,} 404-416.

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{
   
# Create and fit a model using mxExpectationLISREL, and mxFitFunctionML

library(OpenMx)

covData <- matrix(c(0.9223099, 0.1862938, 0.4374359, 0.8959973, 0.9928430, 0.5320662,
                    0.1862938, 0.2889364, 0.3927790, 0.3321639, 0.3371594, 0.4476898,
                    0.4374359, 0.3927790, 1.0069552, 0.6918755, 0.7482155, 0.9013952,
                    0.8959973, 0.3321639, 0.6918755, 1.8059956, 1.6142005, 0.8040448,
                    0.9928430, 0.3371594, 0.7482155, 1.6142005, 1.9223567, 0.8777786,
                    0.5320662, 0.4476898, 0.9013952, 0.8040448, 0.8777786, 1.3997558
                  ), nrow=6, ncol=6, byrow=TRUE, 
                  dimnames=list(paste("v",as.character(1:6),sep=""),paste("v",as.character(1:6),sep="")))

# Create LISREL matrices

mLX <- mxMatrix("Full", values=c(.5, .6, .8, rep(0, 6), .4, .7, .5), 
          name="LX", nrow=6, ncol=2, free=c(TRUE,TRUE,TRUE,rep(FALSE, 6),TRUE,TRUE,TRUE),
          dimnames=list(paste("v",as.character(1:6),sep=""),c("x1","x2")))
mTD <- mxMatrix("Diag", values=c(rep(.2, 6)), 
          name="TD", nrow=6, ncol=6, free=TRUE,
          dimnames=list(paste("v",as.character(1:6),sep=""),paste("v",as.character(1:6),sep="")))
mPH <- mxMatrix("Symm", values=c(1, .3, 1), 
          name="PH", nrow=2, ncol=2, free=c(FALSE, TRUE, FALSE),
          dimnames=list(c("x1","x2"),c("x1","x2")))

# Create a LISREL objective with LX, TD, and PH matrix names

expFunction <- mxExpectationLISREL(LX="LX", TD="TD", PH="PH")

# Create fit function and data

tmpData <- mxData(observed=covData, type="cov", numObs=100)
fitFunction <- mxFitFunctionML()

# Create the model, fit it, and print a summary.

tmpModel <- mxModel(model="exampleModel", mLX, mTD, mPH, expFunction, fitFunction, tmpData)
tmpModelOut <- mxRun(tmpModel)
summary(tmpModelOut)


#--------------------------------------
# Fit factor model with means

require(OpenMx)

data(demoOneFactor)
nvar <- ncol(demoOneFactor)
varnames <- colnames(demoOneFactor)

factorMeans <- mxMatrix("Zero", 1, 1, name="Kappa", dimnames=list("F1", NA))
xIntercepts <- mxMatrix("Zero", nvar, 1, name="TauX", dimnames=list(varnames, NA))
factorLoadings <- mxMatrix("Full", nvar, 1, TRUE, .6, name="LambdaX", labels=paste("lambda", 1:nvar, sep=""), dimnames=list(varnames, "F1"))
factorCovariance <- mxMatrix("Diag", 1, 1, FALSE, 1, name="Phi")
xResidualVariance <- mxMatrix("Diag", nvar, nvar, TRUE, .2, name="ThetaDelta", labels=paste("theta", 1:nvar, sep=""))

liModel <- mxModel(model="LISREL Factor Model",
	factorMeans, xIntercepts, factorLoadings, factorCovariance, xResidualVariance,
	mxExpectationLISREL(LX="LambdaX", PH="Phi", TD="ThetaDelta", TX="TauX", KA="Kappa"),
	mxFitFunctionML(),
	mxData(demoOneFactor, "raw")
)

liRun <- mxRun(liModel)

summary(liRun)


#--------------------------------------
# Fit Joint Ordinal/Continuous Factor Model

require(OpenMx)

# get data (loaded from demo data sets in OpenMx package)
data(jointdata)

# specify ordinal columns as ordered factors
jointdata[,c(2,4,5)] <- mxFactor(jointdata[,c(2,4,5)], 
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

loadings <- mxMatrix("Full", 5, 1,
	free=TRUE, values=1, name="L", dimnames=list(names(jointdata), "Factor1"))

resid <- mxMatrix("Diag", 5, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=.5, name="U")

means <- mxMatrix("Full", 5, 1,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=0, name="M", dimnames=list(names(jointdata), NA))

ident <- mxMatrix("Diag", 1, 1, FALSE, 1, name="I")
zerom <- mxMatrix("Zero", 1, 1, name="Z", dimnames=list("Factor1", NA))

thrFre <- c(TRUE, FALSE, FALSE, rep(TRUE, 5), FALSE)
thrVal <- c(0, NA, NA, -1, 0, 1, -1, 1, NA)
thrLab <- c("z2t1", NA, NA, "z4t1", "z4t2", "z4t3", "z5t1", "z5t2", NA)
thresh <- mxMatrix("Full", 3, 3, free=thrFre, values=thrVal, labels=thrLab, name="T", dimnames=list(c(NA, NA, NA), c("z2", "z4", "z5")))


# run factor model
jointModel1 <- mxModel("ContinuousOrdinalData",
	mxData(jointdata, "raw"),
	loadings, resid, means, ident, zerom, thresh,
	mxFitFunctionML(),
	mxExpectationLISREL(LX="L", TX="M", PH="I", KA="Z", TD="U",
		dimnames=names(jointdata),
		thresholds="T",
		threshnames=c("z2", "z4", "z5"))
)

# Heads up, running this model could take up to 30 seconds.
jointResults1 <- mxRun(jointModel1, suppressWarnings=TRUE)



}
