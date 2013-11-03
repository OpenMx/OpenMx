%
%   Copyright 2007-2013 The OpenMx Project
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

\name{mxPath}
\alias{mxPath}
\alias{MxPath-class}
\alias{print,MxPath-method}
\alias{show,MxPath-method}

\title{Create List of Paths}

\description{
   This function creates a list of paths.
}

\usage{
mxPath(from, to = NA, connect = c("single", "all.pairs", "unique.pairs",
    "all.bivariate", "unique.bivariate"), arrows = 1,
    free = TRUE, values = NA, labels = NA, 
    lbound = NA, ubound = NA, ...)
}

\arguments{
   \item{from}{character vector. These are the sources of the new paths.}
   \item{to}{character vector. These are the sinks of the new paths.}
   \item{connect}{String. Specifies the type of source to sink connection: "single", "all.pairs", "all.bivariate", "unique.pairs", "unique.bivariate". Default value is "single".}
   \item{arrows}{numeric value. Must be either 1 (for single-headed) or 2 (for double-headed arrows).}
   \item{free}{boolean vector. Indicates whether paths are free or fixed.}
   \item{values}{numeric vector. The starting values of the parameters.}
   \item{labels}{character vector. The names of the paths.}
   \item{lbound}{numeric vector. The lower bounds of free parameters.}
   \item{ubound}{numeric vector. The upper bounds of free parameters.}
   \item{...}{Not used.  Allows OpenMx to catch the use of the deprecated \sQuote{all} argument.}
}

\note{
The previous implementation of \sQuote{all} had unsafe features. Its use is now deprecated, and has been replaced by the new mechanism \sQuote{connect} which supports safe and controlled generation of desired combinations of paths. 
}

\details{
The mxPath function creates \link[=MxPath-class]{MxPath} objects. These consist of a list of paths describing the relationships between variables in a model using the RAM modeling approach (McArdle and MacDonald, 1984). Variables are referenced by name, and these names must appear in the \sQuote{manifestVar} and \sQuote{latentVar} arguments of the \link{mxModel} function.

Paths are specified as going "from" one variable (or set of variables) "to" another variable or set of variables using the \sQuote{from} and \sQuote{to} arguments, respectively. If \sQuote{to} is left empty, it will be set to the value of \sQuote{from}.

\sQuote{connect} has five possible connection types: "single", "all.pairs", "all.bivariate", "unique.pairs", "unique.bivariate". The default value is "single". Assuming the values c(\sQuote{a},\sQuote{b},\sQuote{c}) for the \sQuote{to} and \sQuote{from} fields the paths produced by each connection type are as follows:

    \describe{
    \item{"all.pairs":}{(a,a), (a,b), (a,c), (b,a), (b,b), (b,c), (c,a), (c,b), (c,c).} 
    \item{"unique.pairs":}{(a,a), (a,b), (a,c), (b,b), (b,c), (c,c). }
    \item{"all.bivariate":}{(a,b), (a,c), (b,a), (b,c), (c,a), (c,b).} 
    \item{"unique.bivariate":}{(a,b), (a,c), (b,c).} 
    \item{"single":}{(a,a), (b,b), (c,c).}  
    }
    
Multiple variables may be input as a vector of variable names. If the \sQuote{connect} argument is set to "single", then paths are created going from each entry in the \sQuote{from} vector to the corresponding entry in the \sQuote{to} vector. If the \sQuote{to} and \sQuote{from} vectors are of different lengths when the \sQuote{connect} argument is set to "single", the shorter vector is repeated to make the vectors of equal length.

    
The \sQuote{free} argument specifies whether the paths created by the mxPath function are free or fixed parameters. This argument may take either TRUE for free parameters, FALSE for fixed parameters, or a vector of TRUEs and FALSEs to be applied in order to the created paths.

The \sQuote{arrows} argument specifies the type of paths created. A value of 1 indicates a one-headed arrow representing regression. This path represents a regression of the \sQuote{to} variable on the \sQuote{from} variable, such that the arrow points to the \sQuote{to} variable in a path diagram. A value of 2 indicates a two-headed arrow, representing a covariance or variance. If multiple paths are created in the same mxPath function, then the \sQuote{arrows} argument may take a vector of 1s and 2s to be applied to the set of created paths.

The \sQuote{values} is a numeric vectors containing the starting values of the created paths. \sQuote{values} gives a starting value for estimation. The \sQuote{labels} argument specifies the names of the resulting \link[=MxPath-class]{MxPath} object.  The \sQuote{lbound} and \sQuote{ubound} arguments specify lower and upper bounds for the created paths.

}

\value{
Returns a list of paths.
}

\references{
McArdle, J. J. and MacDonald, R. P. (1984). Some algebraic properties of the Reticular Action Model for moment structures. \emph{British Journal of Mathematical and Statistical Psychology, 37,} 234-251.

The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\seealso{
\link{mxMatrix} for a matrix-based approach to path specification; \link{mxModel} for the container in which mxPaths are embedded. More information about the OpenMx package may be found \link[=OpenMx]{here}. 
}


\examples{
# A simple Example: 1 factor Confirmatory Factor Analysis
require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents   <- c("G")
factorModel <- mxModel("One Factor", type="RAM",
      manifestVars = manifests,
      latentVars   = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,free=FALSE, values=1.0),
      mxData(cov(demoOneFactor), type="cov",numObs=500)
)
factorFit <-mxRun(factorModel)
summary(factorFit)

# A more complex example using features of R to compress what would otherwise be a long and error-prone script

myManifest <- sprintf("\%02d", c(1:100))  # list of 100 variable names: "01"  "02"  "03"...
myLatent   <- c("G1", "G2", "G3", "G4", "G5") # the latent variables for the model
# Start building the model: Define its type, and add the manifest and latent variable name lists
model     <- mxModel(type = "RAM", manifestVars = myManifest, latentVars = myLatent)

# Create covariances between the latent variables and add to the model
# Here we use combn to create the covariances
# nb: To create the variances and covariances in one operation you could use 
# expand.grid(myLatent,myLatent) to specify from and to

uniquePairs <- combn(myLatent,2)
covariances <- mxPath(from = uniquePairs[1,], to=uniquePairs[2,], arrows = 2, free = TRUE, values = 1)
model       <- mxModel(model, covariances) 

# Create variances for the latent variables
variances   <- mxPath(from = myLatent, to=myLatent, arrows = 2, free = TRUE, values = 1)
model       <- mxModel(model, variances) # add variances to the model

# Make a list of paths from each packet of 20 manifests to one of the 5 latent variables
# nb: The first loading to each latent is fixed to 1 to scale its variance.
singles <- list()
for (i in 1:5) {
    j <- i*20
    singles <- append(singles, mxPath(
                        from = myLatent[i], to = myManifest[(j - 19):j], 
                        arrows = 1,
                        free   = c(FALSE, rep(TRUE, 19)), 
                        values = c(1, rep(0.75, 19))))
}

model <- mxModel(model, singles) # add single-headed paths to the model

}
