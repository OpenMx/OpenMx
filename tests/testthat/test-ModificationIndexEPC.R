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

# Checks MIs and EPC against a simple model' true value

library(OpenMx)
library(testthat)
context("ModificationIndexEPC")

S <- matrix(c(1, 0.5, 0.5, 1.25), 2, 2, dimnames = list(c("x1", "x2"), c("x1", "x2")))
means <- c(x1 = 0.5, x2 = 1.25)

model <- mxModel("EPC",
    type = "RAM", manifestVars = colnames(S),
    mxPath(
        from = colnames(S), arrows = 2, connect = "unique.pairs", free = FALSE,
        values = S[lower.tri(S, diag = TRUE)]
    ),
    mxPath(from = "one", to = "x1"),
    mxData(S, type = "cov", numObs = 100, means = means)
)
fit <- mxRun(model, silent = TRUE)
mi <- mxMI(fit, matrices = "M")

refit <- mxRun(mi$plusOneParamModels[["EPC.M[1,2]"]], silent = TRUE)
omxCheckCloseEnough(mi$MI.Full[["EPC.M[1,2]"]], fit$output$fit - refit$output$fit, 1e-3)
omxCheckCloseEnough(mi$EPC[["EPC.M[1,2]"]], refit$M$values[1, 2] - fit$M$values[1, 2], 1e-4)
