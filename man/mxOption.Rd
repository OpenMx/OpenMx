\name{mxOption}
\alias{mxOption}

\title{Set or Clear an Optimizer Option}

\description{
The function sets or clears an option that is specific to the optimizer in the back-end.
}

\usage{
mxOption(model, key, value, reset = FALSE)
}

\arguments{
\item{model}{An \link{MxModel} object}
\item{key}{The name of the option.}
\item{value}{The value of the option.}
\item{reset}{If TRUE then reset all options to their defaults.}
}

\details{
Sets an option that is specific to the particular optimizer used in the back-end. 
The name of the option is the \sQuote{key} argument. Use value = NULL to remove 
an existing option. Before the model is submitted to the back-end, all keys and 
values are converted into strings using the \link{as.character} function. To reset
all options to their default values, use reset = TRUE. If reset = TRUE, then 'key'
and 'value' are ignored. Use \code{getOption('mxOptimizerOptions')} to see the default optimizer options.

OpenMx options

\tabular{rccl}{
Calculate Hessian \tab [Yes|No] \tab  \tab calculate the hessian explicitly after optimization.  If \sQuote{No}, a less precise hessian estimate (from the optimizer) is returned. Forced to \sQuote{Yes} if \sQuote{Standard Errors} option is enabled. \cr
Standard Errors \tab [Yes|No] \tab   \tab return standard error estimates from the explicitly calculate hessian.  Forces explicit calculation of the hessian. \cr
Chi-square Confidence Intervals \tab [Yes|No] \tab  \tab calculate confidence intervals based on boundaries of the maximum likelihood space. At present, this assumes a Chi-square distribution for the results of the objective function. \cr
}

NPSOL-specific options

\tabular{rcl}{
Nolist \tab \tab this option suppresses printing of the options \cr
Print level \tab \var{i} \tab the value of \var{i} controls the amount of printout produced by the major iterations \cr
Minor print level \tab \var{i} \tab the value of \var{i} controls the amount of printout produced by the minor iterations \cr
Print file \tab \var{i} \tab for \var{i} > 0 a full log is sent to the file with logical unit number \var{i}. \cr
Summary file \tab \var{i} \tab for \var{i} > 0 a brief log will be output to file \var{i}. \cr
Function Precision \tab \var{r} \tab a measure of accuracy with which \var{f} and \var{c} can be computed. \cr
Infinite Bound Size \tab \var{r} \tab if \var{r} > 0 defines the "infinite" bound bigbnd. \cr 
Feasibility tolerance \tab \var{r} \tab the maximum acceptable absolute violations in linear and nonlinear constraints. \cr
Major iterations \tab \var{i} \tab the maximum number of major iterations before termination. \cr
Verify level \tab [-1:3|Yes|No] \tab see NPSOL manual. \cr
Line search tolerance \tab \var{r} \tab controls the accuracy with which a step is taken. \cr
Derivative Level \tab [0-3] \tab see NPSOL manual. \cr
Hessian \tab [Yes|No] \tab return the transformed Hessian (if \sQuote{No}) or the Hessian itself (if \sQuote{Yes}). \cr
}

}
	
\value{
Returns the model with the optimizer option either set or cleared.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}


\examples{
model <- mxModel()

model <- mxOption(model, "Function Precision", 1e-5)
model <- mxOption(model, "Standard Errors", "Yes")
model <- mxOption(model, "Function Precision", NULL)

}
