library(OpenMx)

symmMatrix <- mxMatrix("Symm", nrow = 3, ncol = 3)
symmMatrix@values[1,2] <- 2
symmMatrix@values[2,1] <- 2 # now required explicitly
omxCheckIdentical(symmMatrix[2,1], symmMatrix[1,2])
