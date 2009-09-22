\name{mxEval}
\alias{mxEval}

\title{Function To Evaluate MxModel Values}

\description{
   This function can be used to evaluate an arbitrary R expression that includes named entities from a \link{MxModel} object, or labels from a \link{MxMatrix} object.
}

\usage{
mxEval(expression, model, compute = FALSE, show = FALSE)
}

\arguments{
   \item{expression}{An arbitrary R expression.}
   \item{model}{The model in which to evaluate the expression.}
   \item{compute}{If TRUE then compute the value of algebra expressions.}
   \item{show}{If TRUE then print the translated expression.}
}

\details{
The argument \sQuote{expression} is an arbitrary R expression.  Any named entities that are used within the R expression are translated into their current value from the model. Any labels from the matrices within the model are translated into their current value from the model. Finally the expression is evaluated and the result is returned.  To enable debugging, the \sQuote{show} argument has been provided.  The most common mistake when using this function is to include named entities in the model that are identical to R function names.  For example, if a model contains a named entity named \sQuote{c}, then the following mxEval call will return an error: \code{mxEval(c(A, B, C), model)}.

If \sQuote{compute} is FALSE, then MxAlgebra expressions returns their current value as they have been computed by the optimization call (using \link{mxRun}). If the \sQuote{compute} argument is TRUE, then MxAlgebra expressions will be calculated in R.  Any references to an objective function that has not yet been calculated will return a 1 x 1 matrix with a value of NA.}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}

\examples{

matrixA <- mxMatrix("Full", nrow = 1, ncol = 1, 
	values = 1, name = "A")
algebraB <- mxAlgebra(A + A, name = "B")

model <- mxModel(matrixA, algebraB)
model <- mxRun(model)
start <- mxEval(-pi * A, model)

\dontrun{
	
mxEval(plot(sin, start, B * pi), model)

# The statement above is equivalent to:

plot(sin, -pi, 2 * pi)

}

}
