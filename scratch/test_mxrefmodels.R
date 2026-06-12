library(OpenMx)

set.seed(42)
x=rnorm(1000, mean=0, sd=1)
y= 0.5*x + rnorm(1000, mean=0, sd=1)
tmpFrame <- data.frame(x, y)
tmpNames <- names(tmpFrame)
wdata <- mxData(tmpFrame, type="raw")

S <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(1,0,0,1),
              free=c(TRUE,FALSE,FALSE,TRUE), labels=c("Vx", NA, NA, "Vy"), name = "S")
A <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(0,1,0,0),
              free=c(FALSE,TRUE,FALSE,FALSE), labels=c(NA, "b", NA, NA), name = "A")
I <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")

expCov <- mxAlgebra(solve(I-A) %*% S %*% t(solve(I-A)), name="expCov")
expFunction <- mxExpectationNormal(covariance="expCov", dimnames=tmpNames)
fitFunction <- mxFitFunctionWLS('DWLS')

tmpModel <- mxModel(model="exampleModel", S, A, I, expCov, expFunction, fitFunction,
                    wdata)

tmpModelOut <- mxRun(tmpModel)
ref <- mxRefModels(tmpModelOut, run=TRUE)

message("--- Summary of the WLS model with reference models ---")
s <- summary(tmpModelOut, refModels=ref)
print(s)
