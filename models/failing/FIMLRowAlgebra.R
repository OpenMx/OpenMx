require(OpenMx)
require(MASS)

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
numRows <- 1000
x <- (mvrnorm(n=numRows, rep(0, 3), expectedcov)) 
dimnames(x) <- list(NULL, c('a','b','c'))

# throw in a few missing values 
#
# x[x>1]<-NA 
# x[4,1]<-NA 
# x

## Now apply the missdmnorm function to the data frame with missing values in it...

rTime <- system.time(inSum <- sum(apply(x, 1, missdmnormIn, mu=expectedmean, sigma=expectedcov), na.rm=TRUE), gcFirst=TRUE)


model <- mxModel()
data <- mxData(data.frame(x), 'raw')

model <- mxModel(model,
    mxMatrix("Symm", 
        values = expectedcov, 
        name = "expectedCov", 
        dimnames = list(c('a','b','c'), c('a','b','c'))),
    mxMatrix("Full",
        name = "dataRow",
        nrow=1, ncol=3, 
        free=FALSE,
        dimnames = list(NULL, c('data.a','data.b','data.c'))),
    mxAlgebra(
        (2*pi)^3 * 1/sqrt(det(expectedCov)) * (t(dataRow) %*% (solve(expectedCov)) %*% dataRow)^(1/2),
        name="rowAlgebra"),
    mxMatrix("Full",
        name="rowResults",
        nrow=numRows, ncol=1,
        values = rep(0, numRows),
        free=FALSE),
    mxAlgebra(sum(rowResults), name="reduceAlgebra"),
    mxRowObjective(rowAlgebra="rowAlgebra", rowResults="rowResults"),
    mxData(data.frame(x), 'raw')
)

# Add the objective function and the data to the model
model <- mxRun(model)

NPSOLOutput <- model@output
outSum <- NPSOLOutput$minimum

omxCheckCloseEnough(inSum, outSum, epsilon = 10 ^ -4)

