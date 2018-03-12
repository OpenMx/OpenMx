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

require(OpenMx)
#mxOption(NULL, "Default optimizer", "NPSOL")

m1 <- mxModel(
  mxMatrix(nrow=1,ncol=1,free=TRUE,name="x", values=.5),
  mxAlgebra(abs(x), name="absX"),
  mxFitFunctionAlgebra("absX"))

m1 <- mxRun(m1)
omxCheckCloseEnough(c(m1$compute$steps[[2]]$output$gradient[1,c('forward','central','backward')]),
                    c(-1, 0, 1), 1e-2)
if (0) {
  # We cannot detect this reliably.
  omxCheckEquals(m1$output$status$code,6)
  omxCheckEquals(summary(m1)$npsolMessage, "The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)")
}

omxCheckError(mxOption(m1, "Gradient algorithm", "forward"), "'Gradient algorithm' is a global option and cannot be set on models.
To change 'Gradient algorithm' globally, use, e.g.:
mxOption(NULL, 'Gradient algorithm', 'forward')")

mxOption(NULL, "Gradient algorithm", "central")
omxCheckEquals(mxComputeGradientDescent()$gradientAlgo, "central")

mxOption(NULL, "Gradient algorithm", "forward")
omxCheckEquals(mxComputeGradientDescent()$gradientAlgo, "forward")
