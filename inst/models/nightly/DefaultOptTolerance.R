#options(error = utils::recover)
library(OpenMx)

numDims <- 2

genOrthogonal<-function(dim) { 
  Q<-MOrthogonal(runif(dim))
  return(Q)
}

# Construct an orthogonal matrix whose first few columns are standardized 'M'
# where columns of 'M' are orthogonal.
# Here "standardized 'M'" means each its columns has length 1.
MOrthogonal<-function(M)
{
  # can set the parameter "tol" of "qr" to decide how small value should be 0
  tmp<-qr(M)
  Q<-qr.Q(tmp,complete=TRUE)
  if(is.vector(M)) { if(Q[1]*M[1]<0) Q<- -Q }
  else { if(Q[1,1]*M[1,1]<0) Q<- - Q }
  return(Q)
}

# adapted from clusterGeneration 1.3.1 by Weiliang Qiu, Harry Joe
genPositiveDefMat <- function(dim, low=-1.4, upp=1.4) {
  u<-matrix(0, dim,dim)
  egvalues <- exp(runif(dim,min=low,max=upp))
  diag(u)<-egvalues #the diagonal elements of u are positive
  Sigma<-u
  if(dim>1)
  { Q<-genOrthogonal(dim) # generate an orthogonal matrix 
    Sigma<-Q%*%u%*%t(Q) # the final positive definite matrix
  }
  Sigma
}

if (0) {
	if (0) {
		mxOption(NULL, "Default optimizer", "NPSOL")
	} else {
		mxOption(NULL, "Default optimizer", "CSOLNP")
	}
}

drift <- 0

for (trial in 1:10) {
  set.seed(trial)
  trueMean <- rnorm(numDims)
  trueCov <- round(genPositiveDefMat(numDims), 6)
  names(trueMean) <- paste("f", 1:numDims, sep="")
  dimnames(trueCov) <- list(paste("f", 1:numDims, sep=""), paste("f", 1:numDims, sep=""))
  
  t1 <- mxModel("test",
                mxData(observed=trueCov, means=trueMean, numObs=100, type="cov"),
                mxMatrix(name="mean", values=rnorm(numDims), nrow=1, ncol=numDims, free=TRUE,
                         dimnames=list(NULL, paste("f", 1:numDims, sep=""))),
                mxMatrix(type="Symm", name="cov", ncol=numDims, nrow=numDims, free=TRUE,
                         labels=paste("v", 1:(numDims*(numDims+1)/2), sep=""),
                         dimnames=dimnames(trueCov)),
                mxExpectationNormal(means="mean", covariance="cov"),
                mxFitFunctionML(),
		mxComputeGradientDescent())
  t1$cov$values <- diag(numDims)
  
  t1Fit <- mxRun(t1, silent=TRUE, suppressWarnings=TRUE)
  
  drift <- drift + sqrt(sum(c(t1Fit$mean$values - trueMean,
			      t1Fit$cov$values - trueCov)^2))
}

print(drift)

omxCheckCloseEnough(drift, 0, 8e-6)
