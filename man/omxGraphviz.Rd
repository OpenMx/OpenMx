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
   Invisibly returns a string containing the model description in graphviz format.
}


\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}
