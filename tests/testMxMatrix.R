library(OpenMx)

symmMatrix <- mxMatrix("Symm", nrow = 3, ncol = 3)
symmMatrix[1,2] <- 2
symmMatrix[2,1] <- 2 # now required explicitly
omxCheckEquals(symmMatrix[2,1], symmMatrix[1,2])
