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


library(OpenMx)

#mxOption(NULL, "Default optimizer", "NPSOL")

# Constrain matrix 'K' to be equal to matrix 'limit'
model <- mxModel(model="con_test",
                 mxMatrix(type="Full", nrow=2, ncol=2, free=TRUE, name="K"),
                 mxMatrix(type="Full", nrow=2, ncol=2, free=FALSE, name="limit", values=1:4),
                 mxConstraint(K == limit, name = "Klimit_equality"),
                 mxAlgebra(min(K), name="minK"),
                 mxFitFunctionAlgebra("minK")
)

fit <- mxRun(model)

omxCheckCloseEnough(fit$output$fit, 1, 1e-4)
omxCheckCloseEnough(c(fit$matrices$K$values), 1:4, 1e-3)
