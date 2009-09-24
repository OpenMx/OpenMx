#
#   Copyright 2007-2009 The OpenMx Project
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
omxCheckEquals(mxEval(D, model), as.matrix(c(1,2,3)))
A <- mxMatrix('Full', 2, 2, values = c(1,2,3,4), byrow = TRUE, name = 'A')
B <- mxAlgebra(A + A, name = 'B')
C <- mxMatrix('Full', 2, 2, labels = c('B[2,2]', 'B[2,1]', 'B[1,2]', 'B[1,1]'), byrow = TRUE, name = 'C')
model <- mxModel('model', A, B, C)
model <- mxRun(model)
omxCheckEquals(mxEval(C, model), matrix(rev(c(1,2,3,4)) * 2, 2, 2, byrow = TRUE))
