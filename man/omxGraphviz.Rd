%
%   Copyright 2007-2021 by the individuals mentioned in the source code history
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

\name{omxGraphviz}
\alias{omxGraphviz}

\title{Show RAM Model in Graphviz Format}

\description{
   The function accepts a RAM style model and outputs a visual representation
   of the model in Graphviz format. The function will output either to a file or 
   to the console. The recommended file extension for an output file is ".dot".
}

\usage{
omxGraphviz(model, dotFilename = "")
}

\arguments{
   \item{model}{An RAM-type model.}
   \item{dotFilename}{The name of the output file.  Use "" to write to console.}
}

\value{
   Invisibly returns a string containing the model description in Graphviz format.
}


\references{
The OpenMx User's guide can be found at \url{https://openmx.ssri.psu.edu/documentation/}.
}
