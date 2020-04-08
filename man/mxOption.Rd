%
%   Copyright 2007-2019 by the individuals mentioned in the source code history
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

\name{mxOption}
\alias{mxOption}

\title{Set or Clear an Optimizer Option}

\description{
The function sets, shows, or clears an option that is specific to the optimizer in the back-end.
}

\usage{
mxOption(model=NULL, key=NULL, value, reset = FALSE)
}

\arguments{
\item{model}{An \link{MxModel} object or NULL}
\item{key}{The name of the option.}
\item{value}{The value of the option.}
\item{reset}{If TRUE then reset all options to their defaults.}
}

\details{
mxOption is used to set, clear, or query an option (given in the \sQuote{key} argument)
in the back-end optimizer. Valid option keys are listed below.

Use value = NULL to remove an existing option. Leaving value blank will return 
the current value of the option specified by \sQuote{key}.

To reset all options to their default values, use \sQuote{reset = TRUE}. 
When reset = TRUE, \sQuote{key} and \sQuote{value} are ignored. 

If the \sQuote{model} argument is set to NULL, the default optimizer option (i.e
those applying to all models by default) will be set.

To see the defaults, use \code{getOption('mxOptions')}.

Before the model is submitted to the back-end, all keys and values are converted into
strings using the \link{as.character} function.

The \dQuote{Default optimizer} option can only be set globally (i.e., with \code{model=NULL}), and not locally (i.e., specifically to a given MxModel).  Although the checkpointing options may be set globally, OpenMx's behavior is only affected by locally set checkpointing options (that is, global checkpointing options are ignored at \link[=mxRun]{runtime}).

Options \dQuote{Gradient step size}, \dQuote{Gradient iterations}, and \dQuote{Function precision} have on-load global defaults of \code{"Auto"}.  If value \code{"Auto"} is in effect for any of these three options at \link[=mxRun]{runtime}, then OpenMx selects a reasonable numerical value in its place.  These automated numerical values are intended to (1) adjust for the limited precision of the algorithm for computing multivariate-normal probability integrals, and (2) calculate accurate numeric derivatives at the optimizer's solution.  If the user replaces \code{"Auto"} with a valid numerical value, then OpenMx uses that value as-is.

\strong{Optimizer specific options}

Some options only affect certain optimizers. Option \dQuote{Analytic Gradients} affects all three optimizers. Option \dQuote{Gradient algorithm} is used by CSOLNP and SLSQP, and ignored by NPSOL. Option \dQuote{Gradient iterations} only affects SLSQP. Option \dQuote{Gradient step size} is used slightly differently by SLSQP and CSOLNP, and is ignored by NPSOL (see \code{\link{mxComputeGradientDescent}()} for details).

If an mxModel contains\link[=MxConstraint]{mxConstraints}, NPSOL is given .4 times the value of the option \dQuote{Feasibility tolerance}. If there are no constraints, NPSOL is given a hard-coded value of 1e-5 (its own native default).

\emph{Note}: Where constraints are present, NPSOL is given 0.4 times the value of the mxOption \dQuote{Feasibility Tolerance}, and this is about a million times bigger than NPSOL's own native default. Values of \dQuote{Feasibility Tolerance} around 1e-5 may be needed to get constraint performance similar to NPSOL's default. Note also that NPSOL's criterion for returning a status code of 0 versus 1 for a given solution depends partly on \dQuote{Optimality tolerance}.

For a block of \code{n} ordinal variables, the maximum number of integration points that OpenMx may use to calculate multivariate-normal probability integrals is given by
\code{mvnMaxPointsA + mvnMaxPointsB*n + mvnMaxPointsC*n*n +
  exp(mvnMaxPointsD + mvnMaxPointsE * n * log(mvnRelEps))}.
Integral approximation is stopped once either \sQuote{mvnAbsEps} or
\sQuote{mvnRelEps} is satisfied.
Use of \sQuote{mvnAbsEps} is deprecated.

The maximum number of major iterations (the option  \dQuote{Major iterations}) 
for optimization for NPSOL  can be specified either by using a 
numeric value (such as 50, 1000, etc) or by specifying a user-defined function.
The user-defined function should accept two arguments as input, the number of 
parameters and the number of constraints, and return a numeric value as output.

OpenMx options

\tabular{rcl}{
Calculate Hessian \tab [Yes | No]       \tab calculate the Hessian explicitly after optimization. \cr
Standard Errors   \tab [Yes | No]       \tab return standard error estimates from the explicitly calculate hessian. \cr
Default optimizer \tab [NPSOL | SLSQP | CSOLNP] \tab the gradient-descent optimizer to use \cr
Number of Threads \tab [0|1|2|...|10|...] \tab number of threads used
for optimization. Default value is taken from the environment variable
OMP_NUM_THREADS or, if that is not set, 2. \cr
Feasibility tolerance \tab \var{r} \tab the maximum acceptable absolute violations in linear and nonlinear constraints. \cr
Optimality tolerance \tab \var{r} \tab the accuracy with which the final iterate approximates a solution to the optimization problem; roughly, the number of reliable significant figures that the fitfunction value should have at the solution. \cr
Gradient algorithm \tab see list \tab finite difference method, either 'forward' or 'central'. \cr
Gradient iterations \tab 1:4 \tab the number of Richardson extrapolation iterations\cr
Gradient step size \tab \var{r} \tab amount of change made to free parameters when numerically calculating gradient\cr
Analytic Gradients \tab [Yes | No] \tab should the optimizer use analytic gradients (if available)?\cr
loglikelihoodScale \tab \var{i} \tab factor by which the loglikelihood is scaled. \cr
Parallel diagnostics \tab [Yes | No] \tab whether to issue diagnostic
messages about use of multiple threads \cr
Nudge zero starts \tab [TRUE | FALSE] \tab Should OpenMx
"nudge" starting values of zero to 0.1 at runtime? \cr
Status OK \tab character vector \tab Status codes that are considered to indicate a successful optimization \cr
Max minutes \tab numeric \tab Maximum backend elapsed time, in minutes
}


NPSOL-specific options

\tabular{rcl}{
Nolist                \tab         \tab this option suppresses printing of the options \cr
Print level           \tab \var{i} \tab the value of \var{i} controls the amount of printout produced by the major iterations \cr
Minor print level     \tab \var{i} \tab the value of \var{i} controls the amount of printout produced by the minor iterations \cr
Print file            \tab \var{i} \tab for \var{i} > 0 a full log is sent to the file with logical unit number \var{i}. \cr
Summary file          \tab \var{i} \tab for \var{i} > 0 a brief log will be output to file \var{i}. \cr
Function precision    \tab \var{r} \tab a measure of accuracy with which the fitfunction and constraint functions can be computed. \cr
Infinite bound size   \tab \var{r} \tab if \var{r} > 0 defines the "infinite" bound bigbnd. \cr 
Major iterations      \tab \var{i} or a function \tab the maximum number of major iterations before termination. \cr
Verify level          \tab [-1:3 | Yes | No]     \tab see NPSOL manual. \cr
Line search tolerance \tab \var{r}    \tab controls the accuracy with which a step is taken. \cr
Derivative level      \tab [0-3]      \tab see NPSOL manual. \cr
Hessian               \tab [Yes | No] \tab return the Hessian (Yes) or the transformed Hessian (No). \cr
Step Limit            \tab \var{r} \tab maximum change in free parameters at first step of linesearch. \cr
}

Checkpointing options

\tabular{rcl}{
Always Checkpoint    \tab [Yes | No]    \tab whether to checkpoint all models during optimization.\cr
Checkpoint Directory \tab    path       \tab the directory into which checkpoint files are written. \cr
Checkpoint Prefix    \tab    string     \tab the string prefix to add to all checkpoint filenames. \cr
Checkpoint Fullpath  \tab    path       \tab overrides the directory and prefix (useful to output to /dev/fd/2) \cr
Checkpoint Units     \tab    see list   \tab the type of units for checkpointing: 'minutes', 'iterations', or 'evaluations'. \cr
Checkpoint Count     \tab    \var{i}    \tab the number of units between checkpoint intervals. \cr
}

Model transformation options

\tabular{rcl}{
Error Checking           \tab [Yes | No] \tab whether model consistency checks are performed in the OpenMx front-end \cr
No Sort Data             \tab            \tab character vector of model names for which FIML data sorting is not performed \cr
RAM Inverse Optimization \tab [Yes | No] \tab whether to enable solve(I - A) optimization \cr
RAM Max Depth            \tab \var{i}    \tab the maximum depth to be used when solve(I - A) optimization is enabled \cr
}

Multivariate normal integration parameters

\tabular{rcl}{
maxOrdinalPerBlock \tab \var{i} \tab maximum number of ordinal variables
to evaluate together \cr
mvnMaxPointsA \tab \var{i} \tab base number of integration points \cr
mvnMaxPointsB \tab \var{i} \tab number of integration points per ordinal variable \cr
mvnMaxPointsC \tab \var{i} \tab number of integration points per squared ordinal variables \cr
mvnMaxPointsD \tab \var{i} \tab see details \cr
mvnMaxPointsE \tab \var{i} \tab see details \cr
mvnAbsEps     \tab \var{i} \tab absolute error tolerance \cr
mvnRelEps     \tab \var{i} \tab relative error tolerance \cr
}

}
    
\value{
If a model is provided, it is returned with the optimizer option either set
or cleared. If value is empty, the current value is returned.
}

\references{
The OpenMx User's guide can be found at http://openmx.ssri.psu.edu/documentation.
}

\seealso{
See \code{\link{mxModel}()}, as almost all uses of \code{mxOption()} are via an mxModel whose options are set or cleared.  See \code{\link{mxComputeGradientDescent}()} for details on how different optimizers are affected by different options. See \link{as.statusCode} for information about the \code{Status OK} option.
}


\examples{
# set the Numbder of Threads (cores to use)
mxOption(key="Number of Threads", value=imxGetNumThreads())

testModel <- mxModel(model = "testModel5") # make a model to use for example
testModel$options   # show the model options (none yet)
options()$mxOptions # list all mxOptions (global settings)

testModel <- mxOption(testModel, "Function precision", 1e-5) # set precision
testModel <- mxOption(testModel, "Function precision", NULL) # clear precision
# N.B. This is model-specific precision (defaults to global setting)

# may optimize for speed
# at cost of not getting standard errors
testModel <- mxOption(testModel, "Calculate Hessian", "No")
testModel <- mxOption(testModel, "Standard Errors"  , "No")

testModel$options # see the list of options you set

}
