require(OpenMx)
require(MASS)

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

model <- mxModel()
model <- mxModel(model, mxMatrix("Symm", expectedcov, name = "covariance"))
model <- mxModel(model, mxMatrix("Zero", name = "means", nrow=1, ncol=3))
data <- mxData(data.frame(x), 'raw')
objective <- mxFIMLObjective(covariance = "covariance", means = "means")

# Add the objective function and the data to the model
model <- mxModel(model, objective, data)

model <- mxRun(model)

NPSOLOutput <- model@output
outSum <- NPSOLOutput$minimum

inSum <- 7.090750e+04 # From SimpleFIML demo
omxCheckCloseEnough(inSum, outSum, epsilon = 10 ^ -4)

