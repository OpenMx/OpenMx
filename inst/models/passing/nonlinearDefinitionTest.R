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

# -----------------

data$z <- mxFactor(z, levels=1:4, labels=c("a",'b','c','d'))
model <- mxModel(model, mxData(data, "raw"))
result <- omxCheckWarning(mxRun(model),
                          "Non Linear Definition.data: definition variable 'z' is a factor; note that it will be treated as integer (as is done by ?unclass). Is this really what you want to do? Really?")
