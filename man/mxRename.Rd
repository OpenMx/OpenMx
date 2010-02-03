\name{mxRename}
\alias{mxRename}

\title{Rename MxModel or a Submodel}

\description{
  This functions renames either the top model or a submodel to a new name.  All internal references to the old model name are replaced with references to the new name.
}

\usage{
mxRename(model, newname, oldname = NA)
}

\arguments{
  \item{model}{a MxModel object.}
  \item{newname}{the new name of the model.}
  \item{oldname}{the name of the target model to rename. If NA then rename top model.}
}


\value{
  Return a \link{mxModel} object with the target model renamed.
}

\references{
The OpenMx User's guide can be found at http://openmx.psyc.virginia.edu/documentation.
}


\examples{

modelA <- mxModel('modelA')
modelB <- mxModel('modelB')
modelC <- mxModel('modelC', modelA, modelB)

# Rename modelC to model1
model1 <- mxRename(modelC, 'model1')

# Rename submodel modelB to model2
model1 <- mxRename(model1, oldname = 'modelB', newname = 'model2')

}
