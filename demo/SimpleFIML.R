require(mnormt)
require(MASS)
require(OpenMx)

# ----------------------------------
# Define the functions.

testFit <- function(objective, startVals=c(), bounds=c(), matList=list(), varList=list(), algList=list(), data=c(), state=c()) {
	return(.Call("callNPSOL", objective, startVals, bounds, matList, varList, algList, data, state));
}

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

x <- (mvrnorm(n=10000, rep(0, 3), expectedcov)) 

# throw in a few missing values 

x[x>1]<-NA 
x[4,1]<-NA 
x

## Now apply the missdmnorm function to the data frame with missing values in it...

#apply(x, 1, missdmnormIn, mu=expectedmean, sigma=expectedcov)

## bingo... -2lnL's

rTime <- system.time(inSum <- sum(apply(x, 1, missdmnormIn, mu=expectedmean, sigma=expectedcov), na.rm=TRUE), gcFirst=TRUE)

model <- mxModel()
model <- mxModel(model, mxMatrix("Symm", expectedcov, name = "covariance"))
model <- mxModel(model, mxMatrix("Zero", name = "means", nrow=1, ncol=3))
data <- mxData(x, 'raw')
objective <- mxFIMLObjective(covariance = "covariance", means = "means")

# Add the objective function and the data to the model
model <- mxModel(model, objective, data)

cTime <- system.time(model <- mxRun(model), gcFirst=TRUE)

NPSOLOutput <- model@output
print(NPSOLOutput)
outSum <- NPSOLOutput$minimum
print(c(inSum, outSum, inSum - outSum, (inSum - outSum) / inSum))
print(c(rTime, cTime))
