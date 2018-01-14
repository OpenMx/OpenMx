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


require(OpenMx)


# ----------------------------------
# Define the functions.

dmnorm <- function(x, mean=rep(0,d), varcov, log=FALSE)
{
  d  <- if(is.matrix(varcov)) ncol(varcov) else 1
  if(d>1 & is.vector(x)) x <- matrix(x, 1, d)
  n  <- if(d==1)  length(x) else nrow(x) 
  X  <- t(matrix(x, nrow=n, ncol=d)) - mean
  Q  <- apply((solve(varcov)%*% X)* X, 2, sum) 
  logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
  logPDF <- as.vector(Q + d*logb(2*pi)+logDet)/(-2)
  if(log) logPDF else exp(logPDF)
}


missdmnormIn <- function(x,mu,sigma) {
   tsel <- !is.na(x)
   if (length(x[tsel]) == 0) return(NA)
   -2 * dmnorm(x[tsel], mu[tsel], sigma[tsel,tsel], log=TRUE)
}

# ----------------------------------
# Test the function with some simulated data.

expectedcov <- rbind(c( 1.0, 0.5, 0.3),
                     c( 0.5, 2.0, 0.2),
                     c( 0.3, 0.2, 3.0))
expectedmean <- c(0,0,0)                     

## simulate some data

x <- (mvtnorm::rmvnorm(n=1000, rep(0, 3), expectedcov)) 
dimnames(x) <- list(NULL, c('a','b','c'))

# throw in a few missing values 

x[x>1]<-NA 
x[4,1]<-NA 
x

## Now apply the missdmnorm function to the data frame with missing values in it...

rTime <- system.time(inSum <- sum(apply(x, 1, missdmnormIn, mu=expectedmean, sigma=expectedcov), na.rm=TRUE), gcFirst=TRUE)


model <- mxModel()
model <- mxModel(model, mxMatrix("Symm", values = expectedcov, name = "covariance",
			dimnames = list(c('a','b','c'), c('a','b','c'))))
model <- mxModel(model, mxMatrix("Zero", name = "means", nrow=1, ncol=3,
			dimnames = list(NULL, c('a','b','c'))))
data <- mxData(data.frame(x), 'raw')
objective <- mxExpectationNormal(covariance = "covariance", means = "means")

# Add the objective function and the data to the model
model <- mxModel(model, objective, data, mxFitFunctionML())

model <- mxRun(model)

NPSOLOutput <- model$output
outSum <- NPSOLOutput$minimum

omxCheckCloseEnough(inSum, outSum, epsilon = 10 ^ -4)

df <- data.frame(foo=1, foo=2, check.names=FALSE)
omxCheckError(mxData(df, 'raw'), "Column names must be unique. Duplicated: 'foo'")
