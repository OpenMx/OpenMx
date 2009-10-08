# SCRIPT: univACE_drop_helper.R
# Timothy  Bates tim.bates@ed.ac.uk
# History:  Thu Oct  8 14:35:02 BST 2009
# OpenMx: http://www.openmx.virginia.com
##########################################

require(OpenMx)

# Add the function library
setParameters <- function(model, labels, free = NA, value = NA) {
   if(missing(model) || !is(model, "MxModel")) {
      stop("The 'model' argument must be a MxModel object")
   }
   if(missing(labels) || typeof(labels) != "character") {
      stop("The 'labels' argument must be a vector of characters")
   }
   if(!is.na(free) && (typeof(free) != "logical" || length(free) != 1)) {
      stop("The 'free' argument must be a single boolean value")
   }
   if(!is.na(value) && (!is.numeric(value) || length(value) != 1)) {
      stop("The 'value' argument must be a single numeric value")
   }
   return(setParametersHelper(model, labels, free, value))
}

setParametersHelper <- function(model, labels, free, value) {
   model@matrices <- lapply(model@matrices, setParametersMatrix, labels, free, value)
   model@submodels <- lapply(model@submodels, setParametersHelper, labels, free, value)
   return(model)
}

setParametersMatrix <- function(matrix, labels, free, value) {
   select <- apply(matrix@labels, c(1,2), function(x) {!is.na(x) && x %in% labels})
   if (!is.na(free)) {
      matrix@free[select] <- free
   }
   if (!is.na(value)) {
      matrix@values[select] <- value
   }
   return(matrix)
}

# Prepare Data
data("twinData", package="OpenMx")
selVars <- c('bmi1','bmi2')
aceVars <- c("A1","C1","E1","A2","C2","E2")
mzfData <- as.matrix(subset(myTwinData, zyg==1, selVars))
dzfData <- as.matrix(subset(myTwinData, zyg==3, selVars))
cov(mzfData, use="pairwise.complete.obs")
cov(dzfData, use="pairwise.complete.obs")

#Fit ACE Model with RawData and Matrices Input
require(OpenMx)
selVars <- c('x','y')
dataMZ <- matrix(c(1,.8,.8,1), nrow = 2, ncol=2, dimnames = list(selVars,selVars))
dataDZ <- matrix(c(1,.5,.5,1), nrow = 2, ncol=2, dimnames = list(selVars,selVars))

modelShare = mxModel("share", 
  mxMatrix("Full",values= .6, free=TRUE, labels='a1', nrow=1, ncol=1, name="a"),
  mxMatrix("Full",values= .6, free=TRUE, labels='c1', nrow=1, ncol=1, name="c"),
  mxMatrix("Full",values= .6, free=TRUE, labels='e1', nrow=1, ncol=1, name="e"),
  mxAlgebra(a * t(a), name="A"),
  mxAlgebra(c * t(c), name="C"),
  mxAlgebra(e * t(e), name="E"),
  mxAlgebra(name="MZcov",
            rbind(cbind(A+C+E, A+C),
                  cbind(A+C,   A+C+E)) 
  ),
  mxAlgebra(name="DZcov",
      rbind(cbind(A+C+E,     .5 %x%A+C),
            cbind(.5 %x%A+C,  A+C+E)) 
  )
)

modelMZ <- mxModel("modelMZ", 
  mxData(observed=dataMZ, type="cov", numObs = 100),
  mxMLObjective(covariance = "share.MZcov",dimnames=selVars)  
)

modelDZ <- mxModel("modelDZ",
  mxData(observed=dataDZ, type="cov", numObs = 100),
  mxMLObjective(covariance = "share.DZcov",dimnames=selVars)
)

model = mxModel("ACE", 
  modelShare, modelMZ, modelDZ,
  mxAlgebra(modelMZ.objective + modelDZ.objective, name="twin"),
  mxAlgebraObjective("twin")
)
fit <- mxRun(model)
summary(fit)

# Now lets drop C using the helper
modelShare = setParameters(modelShare, labels=c("c1"), free = FALSE, value = 0)

# rebuild and run the model
model = mxModel("ACE", 
  modelShare, modelMZ, modelDZ,
  mxAlgebra(modelMZ.objective + modelDZ.objective, name="twin"),
  mxAlgebraObjective("twin")
)
fit <- mxRun(model)
summary(fit)
# hey presto: no c in the fit