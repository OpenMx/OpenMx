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
\item{model}{An \link{MxModel} object or NULL}
\item{key}{The name of the option.}
\item{value}{The value of the option.}
\item{reset}{If TRUE then reset all options to their defaults.}
}

\details{
Sets an option that is specific to the particular optimizer used in the back-end. 
The name of the option is the \sQuote{key} argument. Use value = NULL to remove 
an existing option. Before the model is submitted to the back-end, all keys and 
values are converted into strings using the \link{as.character} function. To reset
all options to their default values, use reset = TRUE. If reset = TRUE, then \sQuote{key}
and \sQuote{value} are ignored. To set the default optimizer options, use the value 
NULL for the \sQuote{model} argument. Use \code{getOption('mxOptions')} to see the default 
optimizer options.

The maximum number of major iterations for the NPSOL optimization (the option 
\dQuote{Major iterations}) can be specified either by using a 
numeric value (such as 50, 1000, etc) or by specifying a user-defined function.
The user-defined function should accept two arguments as input, the number of 
parameters and the number of constraints, and return a numeric value as output.

OpenMx options

\tabular{rcl}{
Calculate Hessian \tab [Yes|No] \tab calculate the hessian explicitly after optimization. \cr
Standard Errors \tab [Yes|No] \tab return standard error estimates from the explicitly calculate hessian. \cr
CI Max Iterations \tab \var{i} \tab the maximum number of retries when calculating confidence intervals. \cr
Default optimizer \tab [NPSOL|CSOLNP] \tab the gradient descent optimizer to use. \cr
}

NPSOL-specific options

\tabular{rcl}{
Nolist \tab \tab this option suppresses printing of the options \cr
Print level \tab \var{i} \tab the value of \var{i} controls the amount of printout produced by the major iterations \cr
Minor print level \tab \var{i} \tab the value of \var{i} controls the amount of printout produced by the minor iterations \cr
Print file \tab \var{i} \tab for \var{i} > 0 a full log is sent to the file with logical unit number \var{i}. \cr
Summary file \tab \var{i} \tab for \var{i} > 0 a brief log will be output to file \var{i}. \cr
Function precision \tab \var{r} \tab a measure of accuracy with which \var{f} and \var{c} can be computed. \cr
Infinite bound size \tab \var{r} \tab if \var{r} > 0 defines the "infinite" bound bigbnd. \cr 
Feasibility tolerance \tab \var{r} \tab the maximum acceptable absolute violations in linear and nonlinear constraints. \cr
Major iterations \tab \var{i} or a function \tab the maximum number of major iterations before termination. \cr
Verify level \tab [-1:3|Yes|No] \tab see NPSOL manual. \cr
Line search tolerance \tab \var{r} \tab controls the accuracy with which a step is taken. \cr
Derivative level \tab [0-3] \tab see NPSOL manual. \cr
Hessian \tab [Yes|No] \tab return the transformed Hessian (if \sQuote{No}) or the Hessian itself (if \sQuote{Yes}). \cr
}

Checkpointing options

\tabular{rl}{
Always Checkpoint \tab "Yes" or "No" whether to checkpoint all models during optimization \cr
Checkpoint Directory \tab the directory where to write checkpoint files \cr
Checkpoint Prefix \tab the string prefix to add to all checkpoint filenames \cr
Checkpoint Units \tab the type of units for checkpointing: 'minutes' or 'iterations' \cr
Checkpoint Count \tab the number of units between checkpoint intervals \cr
Socket Server \tab the server name for sending optimizer state information \cr
Socket Port  \tab the port on the server for sending optimizer state information \cr
Socket Units \tab the type of units: 'minutes' or 'iterations' \cr
Socket Count \tab the number of units between communication to the server \cr
}

Model transformation options

\tabular{rl}{
Error Checking \tab "Yes" or "No" on whether model consistency checks are performed in the OpenMx front-end \cr
No Sort Data   \tab character vector of model names for which FIML data sorting is not performed \cr
RAM Inverse Optimization \tab "Yes" or "No" whether to enable solve(I - A) optimization \cr
RAM Max Depth \tab the maximum depth to be used when solve(I - A) optimization is enabled \cr
}

Multivariate normal integration parameters

\tabular{rl}{
  mvnMaxPointsA \tab base number of integration points \cr
  mvnMaxPointsB \tab number of integration points per row \cr
  mvnMaxPointsC \tab number of integration points per rows^2 \cr
  mvnAbsEps \tab absolute tolerance \cr
  mvnRelEps \tab relative tolerance \cr
}

}
    
\value{
Returns the model with the optimizer option either set or cleared.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxModel} all uses of mxOption are via an mxModel whose options are set or cleared. 
}


\examples{
model <- mxModel() # make a model to use for example
model@options # show the model options (none yet)
options()$mxOptions # list all mxOptions (global settings)

model <- mxOption(model, "Function precision", 1e-5) # set the precision
model <- mxOption(model, "Function precision", NULL) # clear model-specific precision (defaults to global setting)
model <- mxOption(model, "Calculate Hessian", "No")  # may optimize for speed 
model <- mxOption(model, "Standard Errors"  , "No")  # may optimize for speed 
model@options # see the list of options you set

}
