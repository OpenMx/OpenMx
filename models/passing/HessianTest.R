library(OpenMx)
library(numDeriv)
vNames <- c("V1", "V2", "V3", "V4")
dimNames <- list(2)
dimNames[[1]] <- vNames
dimNames[[2]] <- vNames
#dimNames
ObsLD  <- matrix(c( 1,  0,  0,  0,
                   .8,  1,  0,  0,
                   .7, .6,  1,  0,
                   .6, .5, .5,  1), nrow=4, byrow=TRUE)
ObsCov <- ObsLD %*% t(ObsLD)
#ObsCov
dimnames(ObsCov) <- dimNames
model <- mxModel("model",
                 mxMatrix(name="F", type="Full", free=TRUE, nrow=4, ncol=2),
                 mxMatrix(name="U", type="Diag", free=TRUE, nrow=4),
                 mxAlgebra(F %*% t(F) + U, name="PreCov", dimnames <- dimNames),
                 mxData(ObsCov, 'cov', numObs=500),
                 mxMLObjective("PreCov"))
# start values
model@matrices$F@values <- matrix(c(.4, .5, .6, .2, .8, .7, .5, .2), nrow=4, ncol=2)
model@matrices$U@values <- diag(c(.2, .3, .4, .5))                            
# NOTE 10 observed statistics but 12 parameters. model is not identified                               
model <- mxOption(model, "StandardErrors", "Yes")
model <- mxRun(model)                
                                
# now examine the numerically differentiated hessian
fcn <- function(x) {
  F <- matrix(x[1:8], nrow=4, ncol=2)
  U <- diag(x[9:12])
  PreCov <- F %*% t(F) + U
  temp <- ObsCov * solve(PreCov)
  fval <- 499*(log(det(PreCov)) + sum(temp) - 4)
  return(fval)
}
numHess <- hessian(fcn, model@output$estimate, method="Richardson")

omxCheckCloseEnough(numHess, model@output$calculatedHessian, .01)

# use a different set of starting values
model2 <- mxModel("model2",
                 mxMatrix(name="F", type="Full", free=TRUE, nrow=4, ncol=2),
                 mxMatrix(name="U", type="Diag", free=TRUE, nrow=4),
                 mxAlgebra(F %*% t(F) + U, name="PreCov", dimnames <- dimNames),
                 mxData(ObsCov, 'cov', numObs=500),
                 mxMLObjective("PreCov"))
model2@matrices$F@values <- matrix(c(.9, .8, .3, .1, .2, .4, .8, .9), nrow=4, ncol=2)
model2@matrices$U@values <- diag(c(.8, .6, .4, .9))
# again 10 obs stats & 12 parameters
model2 <- mxOption(model2, "StandardErrors", "Yes")
model2 <- mxRun(model2)

# numerical estimate of hessian
numHess2 <- hessian(fcn, model2@output$estimate, method="Richardson")

omxCheckCloseEnough(numHess2, model2@output$calculatedHessian, .01)