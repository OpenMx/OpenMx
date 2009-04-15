\name{mxPath}
\alias{mxPath}
\alias{MxPath}

\title{Function To Create List of Paths}

\description{
   This function creates a list of paths.
}

\usage{
mxPath(from, to = NA, all = FALSE, free = TRUE, 
	arrows = 1, start = NA, end = NA, algebra = NA,
	name = NA, label = NA, boundMax = NA, boundMin = NA,
	ciUpper = NA, ciLower = NA)
}

\arguments{
   \item{from}{character vector. these are the sources of the new paths.}
   \item{to}{character vector. these are the sinks of the new paths.}
   \item{all}{boolean. If TRUE, then connect all sources to all sinks.}
   \item{free}{boolean. Indicates whether paths are free or fixed.}
   \item{arrows}{numeric. Must be either 1 for single-headed or 2 for double-headed arrows.}
   \item{startVal}{numeric. The starting values of the parameters.}
   \item{endVal}{numeric. TODO what does this do?}
   \item{algebra}{character vector. The names of algebras associated with the paths.}
   \item{name}{character vector. The names of the paths.}
   \item{label}{character vector. Descriptive labels for the paths.}
   \item{boundMax}{numeric. TODO what does this do?}
   \item{boundMin}{numeric. TODO what does this do?}
   \item{ciUpper}{numeric. TODO what does this do?}
   \item{ciLower}{numeric. TODO what does this do?}
}

\details{
The mxPath function creates \link{MxPath} objects, which are lists of paths describing the relationships between variables in a model using the RAM modeling approach (McArdle and MacDonald, 1984). Variables are referenced by name, which are included in the 'manifestVar' and 'latentVar' arguments of the \link{mxModel} function. Models using \link{MxPath} objects must be transformed into RAM models using the \link{omxConvertPathModel} function prior to estimation with the \link{mxRAMObjective} function.

Paths are specified as going from one variable or set of variables to another variable or set of variables using the 'from' and 'to' arguments, respectively. Sets of variables may be input as a vector of variable names. If the 'all' argument is set to FALSE, then paths are created going from each entry in the 'from' vector to the corresponding position in the 'to' vector. If the 'to' and 'from' vectors are of different lengths when the 'all' argument is set to FALSE, the shorter vector is repeated to make the vectors of equal length. If the 'all' argument is set to TRUE, all possible paths from the vector of 'from' variables to the vector of 'to' variables are created.

The 'free' argument specifies whether the paths created by the mxPath function are free or fixed parameters. This argument may take either TRUE for free parameters, FALSE for fixed parameters, or a vector of TRUEs and FALSEs to be applied in order to the created paths.

The 'arrows' argument specifies the type of paths created. A value of 1 indicates a one-headed arrow representing regression. This path represents a regression of the 'to' variable on the 'from' variable, such that the arrow points to the 'to' variable in a path diagram. A value of 2 indicates a two-headed arrow, representing a covariance or variance. If multiple paths are created in the same mxPath function, then the 'arrows' argument may take a vector of 1s and 2s to be applied to the set of created paths.

The 'startVal' and 'endVal' are numeric values of the created paths. The 'startVal' gives a starting value for estimation. The 'endVal' may hold a value assigned after estimation. This functionality is not currently used, as \link{MxPath} objects are transformed using the \link{omxConvertPathModel} function prior to estimation, which removes the paths and creates \link{MxMatrix} objects.

The 'algebra' argument specifies optional \link{MxAlgebra} objects associated with the created paths.

The 'name' argument specifies a name of the resulting \link{MxPath} object. This differs from the 'label' argument in that labels refer to the individual paths, while names refer to the \link{MxPath} object as a whole, which may contain multiple paths. 

The 'boundMin' and 'boundMax' arguments specify lower and upper bounds for the created paths.
}

\value{
Returns a list of paths. The set of \link{MxPath} objects included in a model may be transformed into \link{MxMatrix} objects using the \link{omxConvertPathModel} function.
}

\references{
McArdle, J. J. and MacDonald, R. P. (1984). Some algebraic properties of the Reticular Action Model for moment structures. \emph{British Journal of Mathematical and Statistical Psychology, 37,} 234-251.

The Mx manual can be found at http://www.vcu.edu/mx/documentation.html.
}

\examples{

myManifest <- sprintf("\%02d", c(1:100))
myLatent <- c("G1", "G2", "G3", "G4", "G5")
model <- mxModel()

for (i in 1:5) {
    j <- i*20
	singles <- mxPath(from=myLatent[i], to=myManifest[(j-19):j], 
    	arrows=1, free=c(FALSE,rep(TRUE, 19)), 
    	startVal=c(1,rep(0.75,19)))
	model <- mxModel(model, singles)
}

doubles <- mxPath(from=myLatent, all=TRUE, arrows=2,
    free=TRUE, startVal=1)

model <- mxModel(model, doubles)
}
