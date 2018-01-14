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
A <- mxMatrix('Full', 1, 1, values = 1, name = 'A')
B <- mxMatrix('Full', 1, 1, values = 2, name = 'B')
C <- mxMatrix('Full', 1, 1, values = 3, name = 'C')
D <- mxMatrix('Full', 3, 1, labels = c('A[1,1]', 'B[1,1]', 'C[1,1]'), name = 'D')
model <- mxModel('model', A, B, C, D)
model <- mxRun(model)
omxCheckIdentical(mxEval(D, model), as.matrix(c(1,2,3)))
A <- mxMatrix('Full', 2, 2, values = c(1,2,3,4), byrow = TRUE, name = 'A')
B <- mxAlgebra(A + A, name = 'B')
C <- mxMatrix('Full', 2, 2, labels = c('B[2,2]', 'B[2,1]', 'B[1,2]', 'B[1,1]'), byrow = TRUE, name = 'C')
D <- mxMatrix('Full', 1, 2, labels = c('B[2,2]', 'B[2,1]'), byrow = TRUE, name = 'D')
E <- mxAlgebra(t(D), name = 'E')
model <- mxModel('model', A, B, C, D, E)
model <- mxRun(model)
omxCheckIdentical(mxEval(C, model), matrix(rev(c(1,2,3,4)) * 2, 2, 2, byrow = TRUE))
nv <- 2
A <- mxAlgebra(nv, name = 'A')
B <- mxAlgebra(2, name = 'B')
model1 <- mxModel('model1', A)
model2 <- mxModel('model2', B)
model3 <- mxModel('model3', A, B)
model1 <- mxRun(model1)
model2 <- mxRun(model2)
model3 <- mxRun(model3)
omxCheckIdentical(mxEval(A, model1), as.matrix(nv))
omxCheckIdentical(mxEval(B, model2), as.matrix(2))
omxCheckIdentical(mxEval(A, model3), as.matrix(nv))
omxCheckIdentical(mxEval(B, model3), as.matrix(2))
