library(OpenMx)

symmMatrix <- new("SymmMatrix", nrow=3, ncol=3)
symmMatrix[1,2] <- 2
omxCheckEquals(symmMatrix[2,1], symmMatrix[1,2])
