#preamble of some kind

library(OpenMx)

#sample size
n <- 500

set.seed(10)

# generate data
x <- rnorm(n, 3, 2)
z <- rep(1:4, each=n/4)
y <- (.8^z) * x + rnorm(n)

data <- data.frame(x,y,z)

model <- mxModel("Non Linear Definition",
  mxData(data, "raw"),
  mxMatrix("Diag", 2, 2, T, 1, name="S"),
  mxMatrix("Full", 2, 2, 
    free=c(F, F, F, F), 
    values=c(1, 0, 0, 1),  
    labels=c(NA, "bexp[1,1]", NA, NA),
    name="IA"),
  mxMatrix("Full", 1, 1, T, 0.7, "beta1", name="B"),
  mxMatrix("Full", 1, 1, F, 0, "data.z", name="D"),
  mxMatrix("Full", 1, 2, T, 0, c("mu_x", "beta0"), name="M"),
  mxAlgebra(B ^ D, name="bexp"),
  mxAlgebra(M %*% t(IA), name="mu"),
  mxAlgebra(IA %*% S %*% t(IA), name="sigma"),
  mxExpectationNormal("sigma", "mu", dimnames=c("x", "y")),
  mxFitFunctionML()
)

results <- mxRun(model)

summary(results)

check <- nls(y ~ b0 + (b1 ^ z) * x, start=list(b0=0, b1=0.7))

#beta0
omxCheckCloseEnough(results$output$estimate[5], 
  summary(check)$parameters[1], 0.01)

#beta1
omxCheckCloseEnough(results$output$estimate[3], 
  summary(check)$parameters[2], 0.01)
