\name{mxEvaluate}
\alias{mxEvaluate}

\title{Function To Evaluate MxModel Values}

\description{
   This function can be used to evaluate an arbitrary R expression that includes named entities from a \link{MxModel} object.
}

\usage{
mxEvaluate(expression, model, show = FALSE)
}

\arguments{
   \item{expression}{An arbitrary R expression.}
   \item{model}{The model in which to evaluate the expression.}
   \item{show}{If TRUE then print the translated expression.}
}

\details{
The argument \sQuote{expression} is an arbitrary R expression.  Any named entities that are used within the R expression are translated into their current value from the model. Finally the expression is evaluated and the result is returned.  To enable debugging, the \sQuote{show} argument has been provided.  The most common mistake when using this function is to include named entities in the model that are identical to R function names.  For example, if a model contains a named entity named \sQuote{c}, then the following mxEvaluate call will return an error: \code{mxEvaluate(c(A, B, C), model)}.}

\examples{

matrixA <- mxMatrix("Full", nrow = 1, ncol = 1, 
	values = "1", name = "A")
algebraB <- mxAlgebra(A + A, name = "B")

model <- mxModel(matrixA, algebraB)
model <- mxRun(model)
start <- mxEvaluate(-pi * A, model)

\dontrun{
	
mxEvaluate(plot(sin, start, B * pi), model)

# The statement above is equivalent to:

plot(sin, -pi, 2 * pi)

}

}
