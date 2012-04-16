%
%   Copyright 2007-2012 The OpenMx Project
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

\name{mxCompare}
\alias{mxCompare}

\title{Example 500-row dataset with 12 generated variables}

\description{
12 columns of generated data: x1 x2 x3 x4 x5 x6 y1 y2 y3 z1 z2 z3 – each with 500 rows.

The x variables intercorrelate around .6 with each other.

The y variables intercorrelate around .5 with each other, and correlate around .3 with the X vars.

There are three ordinal variables, z1, z2, and z3.

The data are used in some OpenMx examples, especially confirmatory factor analysis.

There are no missing data.
}

\usage{
data("myFADataRaw")
}

\seealso{
# \code{\link{mxModel}};  \code{\link{options}} (use options('mxOptions') to see all the OpenMx-specific options) 
}

\examples{
# load the datafile
data(myFADataRaw)
# extract just the x variables
manifests <- paste("x", 1:6, sep="")
df= myFADataRaw[,manifests]
# Build a 1-factor model of these data
latents <- c("G1")
fit1 <- mxModel("One Factor", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from = latents, to=manifests, labels=paste("path_G1_",manifests, sep="")),
      # residual variance
	  mxPath(from = manifests, arrows = 2, values=.35, labels=paste("residual_",manifests, sep="")), 
      # Scale latent variable to 1
	  mxPath(from = latents, arrows = 2, free = FALSE, values = 1),
	  # Include means for raw data
	  mxPath(from = "one", to = manifests, values = colMeans(df), labels=paste("mean_", manifests, sep="")),
      mxData(df, type = "raw")
)
fit1 <- mxRun(fit1);
# See what thou hath wrought…
summary(fit1)
# omxGraphviz(model=fit1)
}
