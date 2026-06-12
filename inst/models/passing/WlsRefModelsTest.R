#
#   Copyright 2007-2026 by the individuals mentioned in the source code history
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

#------------------------------------------------------------------------------
# Author: Antigravity AI (Pair Programming)
# Date: 2026-06-12
# Filename: WlsRefModelsTest.R
# Purpose: Test WLS reference models (mxRefModels) and summary fit indices (CFI/TLI/RMSEA)
#------------------------------------------------------------------------------

library(OpenMx)

# 1. Generate data
set.seed(42)
x <- rnorm(1000, mean=0, sd=1)
y <- 0.5 * x + rnorm(1000, mean=0, sd=1)
tmpFrame <- data.frame(x, y)
tmpNames <- names(tmpFrame)
wdata <- mxData(tmpFrame, type="raw")

# 2. Define a target model with df > 0 (by fixing path b to 0, which makes it fit like the independence model)
S <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(1,0,0,1),
              free=c(TRUE,FALSE,FALSE,TRUE), labels=c("Vx", NA, NA, "Vy"), name = "S")
A <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(0,0,0,0),
              free=c(FALSE,FALSE,FALSE,FALSE), labels=c(NA, "b", NA, NA), name = "A")
I <- mxMatrix(type="Iden", nrow=2, ncol=2, name="I")

expCov <- mxAlgebra(solve(I-A) %*% S %*% t(solve(I-A)), name="expCov")
expFunction <- mxExpectationNormal(covariance="expCov", dimnames=tmpNames)
fitFunction <- mxFitFunctionWLS("DWLS")

tmpModel <- mxModel(model="exampleModel", S, A, I, expCov, expFunction, fitFunction, wdata)
tmpModelOut <- mxRun(tmpModel)

# 3. Create reference models
ref <- mxRefModels(tmpModelOut, run=TRUE)

# Assertions on reference models
omxCheckEquals(names(ref), c("Saturated", "Independence"))
omxCheckEquals(class(ref$Saturated$fitfunction)[1], "MxFitFunctionWLS")
omxCheckEquals(class(ref$Independence$fitfunction)[1], "MxFitFunctionWLS")

# Saturated model fit for continuous WLS should be 0
omxCheckCloseEnough(ref$Saturated$output$chi, 0, 1e-5)
# Independence model fit should match target model fit because the target model had the path fixed to 0
omxCheckCloseEnough(ref$Independence$output$chi, tmpModelOut$output$chi, 1e-5)

# 4. Generate summary
s <- summary(tmpModelOut, refModels=ref)

# Assertions on summary fit metrics
omxCheckCloseEnough(s$Minus2LogLikelihood, tmpModelOut$output$chi, 1e-5)
omxCheckCloseEnough(s$SaturatedLikelihood, 0, 1e-5)
omxCheckCloseEnough(s$IndependenceLikelihood, ref$Independence$output$chi, 1e-5)

# Fit indices: since the target model has the path fixed to 0, it behaves exactly like the independence model.
# Therefore, CFI and TLI should be 0.
omxCheckCloseEnough(s$CFI, 0, 1e-5)
omxCheckCloseEnough(s$TLI, 0, 1e-5)

# RMSEA should be valid/non-NA
omxCheckTrue(!is.na(s$RMSEA))
omxCheckCloseEnough(s$RMSEA, 0.3894752, 1e-4)

# 5. Define a target model with df = 0 (fully saturated)
A_sat <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values=c(0,0.5,0,0),
              free=c(FALSE,TRUE,FALSE,FALSE), labels=c(NA, "b", NA, NA), name = "A")
tmpModelSat <- mxModel(model="satModel", S, A_sat, I, expCov, expFunction, fitFunction, wdata)
tmpModelSatOut <- mxRun(tmpModelSat)
refSat <- mxRefModels(tmpModelSatOut, run=TRUE)
s_sat <- summary(tmpModelSatOut, refModels=refSat)

# CFI should be 1 for a saturated/just-identified model
omxCheckCloseEnough(s_sat$CFI, 1, 1e-5)
# TLI and RMSEA should be NA (since df = 0)
omxCheckTrue(is.na(s_sat$TLI))
omxCheckTrue(is.na(s_sat$RMSEA))
