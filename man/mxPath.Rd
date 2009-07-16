\name{mxPath}
\alias{mxPath}
\alias{MxPath}

\title{Function To Create List of Paths}

\description{
   This function creates a list of paths.
}

\usage{
mxPath(from, to = NA, all = FALSE, arrows = 1,
	values = NA, free = TRUE, labels = NA, 
	lbound = NA, ubound = NA)
}

\arguments{
   \item{from}{character vector. these are the sources of the new paths.}
   \item{to}{character vector. these are the sinks of the new paths.}
   \item{all}{boolean value. If TRUE, then connect all sources to all sinks.}
   \item{arrows}{numeric value. Must be either 1 for single-headed or 2 for double-headed arrows.}
   \item{values}{numeric vector. The starting values of the parameters.}
   \item{free}{boolean vector. Indicates whether paths are free or fixed.}
   \item{labels}{character vector. The names of the paths.}
   \item{lbound}{numeric vector. The lower bounds of free parameters.}
   \item{ubound}{numeric vector. The upper bounds of free parameters.}   
}

\details{
The mxPath function creates \link{MxPath} objects, which are lists of paths describing the relationships between variables in a model using the RAM modeling approach (McArdle and MacDonald, 1984). Variables are referenced by name, which are included in the 'manifestVar' and 'latentVar' arguments of the \link{mxModel} function.

Paths are specified as going from one variable or set of variables to another variable or set of variables using the 'from' and 'to' arguments, respectively. Sets of variables may be input as a vector of variable names. If the 'all' argument is set to FALSE, then paths are created going from each entry in the 'from' vector to the corresponding position in the 'to' vector. If the 'to' and 'from' vectors are of different lengths when the 'all' argument is set to FALSE, the shorter vector is repeated to make the vectors of equal length. If the 'all' argument is set to TRUE, all possible paths from the vector of 'from' variables to the vector of 'to' variables are created.

The 'free' argument specifies whether the paths created by the mxPath function are free or fixed parameters. This argument may take either TRUE for free parameters, FALSE for fixed parameters, or a vector of TRUEs and FALSEs to be applied in order to the created paths.

The 'arrows' argument specifies the type of paths created. A value of 1 indicates a one-headed arrow representing regression. This path represents a regression of the 'to' variable on the 'from' variable, such that the arrow points to the 'to' variable in a path diagram. A value of 2 indicates a two-headed arrow, representing a covariance or variance. If multiple paths are created in the same mxPath function, then the 'arrows' argument may take a vector of 1s and 2s to be applied to the set of created paths.

The 'values' is a numeric vectors containing the starting values of the created paths. 'values' gives a starting value for estimation. The 'labels' argument specifies the names of the resulting \link{MxPath} object.  The 'lbound' and 'ubound' arguments specify lower and upper bounds for the created paths.
}

\value{
Returns a list of paths.
}

\references{
McArdle, J. J. and MacDonald, R. P. (1984). Some algebraic properties of the Reticular Action Model for moment structures. \emph{British Journal of Mathematical and Statistical Psychology, 37,} 234-251.

The Mx manual can be found at http://www.vcu.edu/mx/documentation.html.
}

\examples{

myManifest <- sprintf("\%02d", c(1:100))
myLatent <- c("G1", "G2", "G3", "G4", "G5")
model <- mxModel(type = "RAM", manifestVars = myManifest, latentVars = myLatent)

singles <- list()
for (i in 1:5) {
    j <- i*20
    singles <- c(singles, mxPath(from = myLatent[i], 
                            to = myManifest[(j - 19) : j], 
                            arrows = 1,
                            free = c(FALSE, rep(TRUE, 19)), 
                            values = c(1, rep(0.75, 19))))
}
model <- mxModel(model, singles)

doubles <- mxPath(from = myLatent, all = TRUE, arrows = 2,
    free = TRUE, values = 1)

model <- mxModel(model, doubles)
}
