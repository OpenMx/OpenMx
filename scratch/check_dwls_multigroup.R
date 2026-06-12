library(OpenMx)

# Simulate some data
set.seed(42)
x=rnorm(1000, mean=0, sd=1)
y= 0.5*x + rnorm(1000, mean=0, sd=1)
tmpFrame <- data.frame(x, y)
tmpNames <- names(tmpFrame)
wdata <- mxData(tmpFrame, type="raw")

# Define the matrices
S <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(1,0,0,1),
              free=c(TRUE,FALSE,FALSE,TRUE), labels=c("Vx", NA, NA, "Vy"), name = "S")
A <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(0,1,0,0),
              free=c(FALSE,TRUE,FALSE,FALSE), labels=c(NA, "b", NA, NA), name = "A")
I <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")

# Define the expectation
expCov <- mxAlgebra(solve(I-A) %*% S %*% t(solve(I-A)), name="expCov")
expFunction <- mxExpectationNormal(covariance="expCov", dimnames=tmpNames)

# Choose a fit function (DWLS)
fitFunction <- mxFitFunctionWLS('DWLS')

# Define the model
tmpModel <- mxModel(model="exampleModel", S, A, I, expCov, expFunction, fitFunction,
                    wdata, mxCI("A"))

tmpModel2 <- mxModel(tmpModel, name="tmp2")
twoGroup <- mxModel("two", tmpModel, tmpModel2, mxFitFunctionMultigroup(c("exampleModel","tmp2")))
twoGroup <- mxRun(twoGroup)

message("--- DWLS Multigroup Model output values ---")
cat("chi: ", twoGroup$output$chi, "\n")
cat("chiDoF: ", twoGroup$output$chiDoF, "\n")
cat("chiM: ", twoGroup$output$chiM, "\n")
cat("chiMV: ", twoGroup$output$chiMV, "\n")
cat("chiMadjust: ", twoGroup$output$chiMadjust, "\n")
cat("chiMVadjust: ", twoGroup$output$chiMVadjust, "\n")
cat("chiDoFstar: ", twoGroup$output$chiDoFstar, "\n")
cat("fitUnits: ", twoGroup$output$fitUnits, "\n")
