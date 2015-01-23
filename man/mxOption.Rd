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

\name{mxOption}
\alias{mxOption}

\title{Set or Clear an Optimizer Option}

\description{
The function sets, shows, or clears an option that is specific to the optimizer in the back-end.
}

\usage{
mxOption(model, key, value, reset = FALSE)
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


The maximum number of major iterations (the option  \dQuote{Major iterations}) 
for optimization for NPSOL  can be specified either by using a 
numeric value (such as 50, 1000, etc) or by specifying a user-defined function.
The user-defined function should accept two arguments as input, the number of 
parameters and the number of constraints, and return a numeric value as output.

OpenMx options

\tabular{rcl}{
Number of Threads \tab \var{i}          \tab the number of processor cores to use. Use \code{\link[parallel]{detectCores}}()
 to find how many are available. \cr
Calculate Hessian \tab [Yes | No]       \tab calculate the Hessian explicitly after optimization. \cr
Standard Errors   \tab [Yes | No]       \tab return standard error estimates from the explicitly calculate hessian. \cr
CI Max Iterations \tab \var{i}          \tab the maximum number of retries when calculating confidence intervals. \cr
Default optimizer \tab [NPSOL | CSOLNP] \tab the gradient descent optimizer to use \cr
Number of Threads \tab [0|1|2|...|10|...] \tab number of threads used for optimization.  This is how parallelism works. Default value of 0 uses \code{detectCores() - 1}. \cr
}

CSOLNP is our new, experimental optimizer. There are known bugs. See the
release note for details.

NPSOL-specific options

\tabular{rcl}{
Nolist                \tab         \tab this option suppresses printing of the options \cr
Print level           \tab \var{i} \tab the value of \var{i} controls the amount of printout produced by the major iterations \cr
Minor print level     \tab \var{i} \tab the value of \var{i} controls the amount of printout produced by the minor iterations \cr
Print file            \tab \var{i} \tab for \var{i} > 0 a full log is sent to the file with logical unit number \var{i}. \cr
Summary file          \tab \var{i} \tab for \var{i} > 0 a brief log will be output to file \var{i}. \cr
Function precision    \tab \var{r} \tab a measure of accuracy with which \var{f} and \var{c} can be computed. \cr
Infinite bound size   \tab \var{r} \tab if \var{r} > 0 defines the "infinite" bound bigbnd. \cr 
Feasibility tolerance \tab \var{r} \tab the maximum acceptable absolute violations in linear and nonlinear constraints. \cr
Major iterations      \tab \var{i} or a function \tab the maximum number of major iterations before termination. \cr
Verify level          \tab [-1:3 | Yes | No]     \tab see NPSOL manual. \cr
Line search tolerance \tab \var{r}    \tab controls the accuracy with which a step is taken. \cr
Derivative level      \tab [0-3]      \tab see NPSOL manual. \cr
Hessian               \tab [Yes | No] \tab return the Hessian (Yes) or the transformed Hessian (No). \cr
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
mvnMaxPointsA \tab \var{i} \tab base number of integration points \cr
mvnMaxPointsB \tab \var{i} \tab number of integration points per row \cr
mvnMaxPointsC \tab \var{i} \tab number of integration points per rows^2 \cr
mvnAbsEps     \tab \var{i} \tab absolute tolerance \cr
mvnRelEps     \tab \var{i} \tab relative tolerance \cr
}

}
    
\value{
If a model is provided, it is returned with the optimizer option either set
or cleared. If value is empty, the current value is returned.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxModel} all uses of mxOption are via an mxModel whose options are set or cleared. 
}


\examples{
# set the Numbder of Threads (cores to use)
mxOption(NULL, "Number of Threads", detectCores() - 1)

testModel <- mxModel(model = "testModel") # make a model to use for example
testModel$options   # show the model options (none yet)
options()$mxOptions # list all mxOptions (global settings)

testModel <- mxOption(testModel, "Function precision", 1e-5) # set the precision
testModel <- mxOption(testModel, "Function precision", NULL) # clear model-specific precision (defaults to global setting)
testModel <- mxOption(testModel, "Calculate Hessian", "No")  # may optimize for speed 
testModel <- mxOption(testModel, "Standard Errors"  , "No")  # may optimize for speed 
testModel$options # see the list of options you set

}
