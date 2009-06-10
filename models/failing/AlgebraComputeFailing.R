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

A <- mxMatrix(values = runif(25), nrow = 5, ncol = 5, name = 'A')
B <- mxMatrix(values = runif(25), nrow = 5, ncol = 5, name = 'B')

model <- mxModel(A, B)

# Insert failing tests
model <- mxModel(model, mxAlgebra(A ^ B, name = 'test3'))
model <- mxModel(model, mxAlgebra(A %&% B, name = 'test7'))
model <- mxModel(model, mxAlgebra(cbind(A,B), name = 'test12'))
model <- mxModel(model, mxAlgebra(rbind(A,B), name = 'test13'))
model <- mxModel(model, mxAlgebra(det(A), name = 'test14'))

model <- mxRun(model)

# Check failing tests
omxCheckCloseEnough(model[['test3']]@result, A@values ^ B@values, 0.001)
omxCheckCloseEnough(model[['test7']]@result, A@values %&% B@values, 0.001)
omxCheckCloseEnough(model[['test12']]@result, cbind(A@values, B@values), 0.001)
omxCheckCloseEnough(model[['test13']]@result, rbind(A@values, B@values), 0.001)
omxCheckCloseEnough(model[['test14']]@result, det(A@values), 0.001)

