#
#   Copyright 2007-2018 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


#options(error = utils::recover)
#library(mvtnorm)
library(OpenMx)

#set.seed(1)

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

R.dmvnorm <- function (x, mean, sigma, log=FALSE)
{
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    distval <- mahalanobis(x, center = mean, cov = sigma)
    logdet <- sum(log(eigen(sigma, symmetric=TRUE,
                                   only.values=TRUE)$values))
    logretval <- -(ncol(x)*log(2*pi) + logdet + distval)/2
    if(log) return(logretval)
    exp(logretval)
}

dmvnorm <- function (x, mean, sigma, log=FALSE)
{
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                     check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    logretval <- apply(x, 1, function(loc) imxDmvnorm(loc, mean, sigma))
    if(log) return(logretval)
    exp(logretval)
}

for (dim in 1:3) {
  sigma <- genPositiveDefMat(dim)
  center <- rnorm(dim)
  loc <- rnorm(dim)

  stopifnot(all.equal(dmvnorm(loc, center, sigma, log=TRUE),
                      R.dmvnorm(loc, center, sigma, log=TRUE)))
}

sigma <- genPositiveDefMat(4)
center <- rnorm(4)
loc <- matrix(rnorm(4*8), nrow=8)
stopifnot(all.equal(R.dmvnorm(loc, center, sigma),
                    dmvnorm(loc, center, sigma)))
