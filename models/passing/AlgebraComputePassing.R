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

# Insert passing tests
model <- mxModel(model, mxAlgebra(A, name = 'test0'))
model <- mxModel(model, mxAlgebra(solve(A), name = 'test1'))
model <- mxModel(model, mxAlgebra(t(A), name = 'test2'))
model <- mxModel(model, mxAlgebra(A ^ B, name = 'test3'))
model <- mxModel(model, mxAlgebra(t(A) ^ B, name = 'test3a'))
model <- mxModel(model, mxAlgebra(A ^ t(B), name = 'test3b'))
model <- mxModel(model, mxAlgebra(A %*% B, name = 'test4'))
model <- mxModel(model, mxAlgebra(t(A) %*% B, name = 'test4a'))
model <- mxModel(model, mxAlgebra(A %*% t(B), name = 'test4b'))
model <- mxModel(model, mxAlgebra(A * B, name = 'test5'))
model <- mxModel(model, mxAlgebra(t(A) * B, name = 'test5a'))
model <- mxModel(model, mxAlgebra(A * t(B), name = 'test5b'))
model <- mxModel(model, mxAlgebra(A %x% B, name = 'test6'))
model <- mxModel(model, mxAlgebra(t(A) %x% B, name = 'test6a'))
model <- mxModel(model, mxAlgebra(A %x% t(B), name = 'test6b'))
model <- mxModel(model, mxAlgebra(A %&% B, name = 'test7'))
model <- mxModel(model, mxAlgebra(A / B, name = 'test8'))
model <- mxModel(model, mxAlgebra(t(A) / B, name = 'test8a'))
model <- mxModel(model, mxAlgebra(A / t(B), name = 'test8b'))
model <- mxModel(model, mxAlgebra(A + B, name = 'test9'))
model <- mxModel(model, mxAlgebra(t(A) + B, name = 'test9a'))
model <- mxModel(model, mxAlgebra(A + t(B), name = 'test9b'))
model <- mxModel(model, mxAlgebra(A - B, name = 'test10'))
model <- mxModel(model, mxAlgebra(t(A) - B, name = 'test10a'))
model <- mxModel(model, mxAlgebra(A - t(B), name = 'test10b'))
model <- mxModel(model, mxAlgebra(-A, name = 'test11'))
model <- mxModel(model, mxAlgebra(-t(A), name = 'test11a'))
model <- mxModel(model, mxAlgebra(t(-A), name = 'test11b'))
model <- mxModel(model, mxAlgebra(cbind(A,B), name = 'test12'))
model <- mxModel(model, mxAlgebra(cbind(t(A),B), name = 'test12a'))
model <- mxModel(model, mxAlgebra(cbind(A,t(B)), name = 'test12b'))
model <- mxModel(model, mxAlgebra(rbind(A,B), name = 'test13'))
model <- mxModel(model, mxAlgebra(rbind(t(A),B), name = 'test13a'))
model <- mxModel(model, mxAlgebra(rbind(A,t(B)), name = 'test13b'))
# test14 is failing
model <- mxModel(model, mxAlgebra(tr(A), name = 'test15'))
model <- mxModel(model, mxAlgebra(sum(A,B), name = 'test16'))
model <- mxModel(model, mxAlgebra(prod(A,B), name = 'test17'))
model <- mxModel(model, mxAlgebra(max(A,B), name = 'test18'))
model <- mxModel(model, mxAlgebra(min(A,B), name = 'test19'))
model <- mxModel(model, mxAlgebra(abs(A), name = 'test20'))
model <- mxModel(model, mxAlgebra(cos(A), name = 'test21'))
model <- mxModel(model, mxAlgebra(cosh(A), name = 'test22'))
model <- mxModel(model, mxAlgebra(sin(A), name = 'test23'))
model <- mxModel(model, mxAlgebra(sinh(A), name = 'test24'))
model <- mxModel(model, mxAlgebra(tan(A), name = 'test25'))
model <- mxModel(model, mxAlgebra(tanh(A), name = 'test26'))
model <- mxModel(model, mxAlgebra(exp(A), name = 'test27'))
model <- mxModel(model, mxAlgebra(log(A), name = 'test28'))
model <- mxModel(model, mxAlgebra(sqrt(A), name = 'test29'))

model <- mxRun(model)

# Check passing tests
omxCheckCloseEnough(model[['test0']]@result, A@values, 0.001)
omxCheckCloseEnough(model[['test1']]@result, solve(A@values), 0.001)
omxCheckCloseEnough(model[['test2']]@result, t(A@values), 0.001)
omxCheckCloseEnough(model[['test3']]@result, A@values ^ B@values, 0.001)
omxCheckCloseEnough(model[['test3a']]@result, t(A@values) ^ B@values, 0.001)
omxCheckCloseEnough(model[['test3b']]@result, A@values ^ t(B@values), 0.001)
omxCheckCloseEnough(model[['test4']]@result, A@values %*% B@values, 0.001)
omxCheckCloseEnough(model[['test4a']]@result, t(A@values) %*% B@values, 0.001)
omxCheckCloseEnough(model[['test4b']]@result, A@values %*% t(B@values), 0.001)
omxCheckCloseEnough(model[['test5']]@result, A@values * B@values, 0.001)
omxCheckCloseEnough(model[['test5a']]@result, t(A@values) * B@values, 0.001)
omxCheckCloseEnough(model[['test5b']]@result, A@values * t(B@values), 0.001)
omxCheckCloseEnough(model[['test6']]@result, A@values %x% B@values, 0.001)
omxCheckCloseEnough(model[['test6a']]@result, t(A@values) %x% B@values, 0.001)
omxCheckCloseEnough(model[['test6b']]@result, A@values %x% t(B@values), 0.001)
omxCheckCloseEnough(model[['test7']]@result, A@values %&% B@values, 0.001)
omxCheckCloseEnough(model[['test8']]@result, A@values / B@values, 0.001)
omxCheckCloseEnough(model[['test8a']]@result, t(A@values) / B@values, 0.001)
omxCheckCloseEnough(model[['test8b']]@result, A@values / t(B@values), 0.001)
omxCheckCloseEnough(model[['test9']]@result, A@values + B@values, 0.001)
omxCheckCloseEnough(model[['test9a']]@result, t(A@values) + B@values, 0.001)
omxCheckCloseEnough(model[['test9b']]@result, A@values + t(B@values), 0.001)
omxCheckCloseEnough(model[['test10']]@result, A@values - B@values, 0.001)
omxCheckCloseEnough(model[['test10a']]@result, t(A@values) - B@values, 0.001)
omxCheckCloseEnough(model[['test10b']]@result, A@values - t(B@values), 0.001)
omxCheckCloseEnough(model[['test11']]@result, -A@values, 0.001)
omxCheckCloseEnough(model[['test11a']]@result, -t(A@values), 0.001)
omxCheckCloseEnough(model[['test11b']]@result, t(-A@values), 0.001)
omxCheckCloseEnough(model[['test12']]@result, cbind(A@values, B@values), 0.001)
omxCheckCloseEnough(model[['test12a']]@result, cbind(t(A@values), B@values), 0.001)
omxCheckCloseEnough(model[['test12b']]@result, cbind(A@values, t(B@values)), 0.001)
omxCheckCloseEnough(model[['test13']]@result, rbind(A@values, B@values), 0.001)
omxCheckCloseEnough(model[['test13a']]@result, rbind(t(A@values), B@values), 0.001)
omxCheckCloseEnough(model[['test13b']]@result, rbind(A@values, t(B@values)), 0.001)
# test14 is failing
omxCheckCloseEnough(model[['test15']]@result, tr(A@values), 0.001)
omxCheckCloseEnough(model[['test16']]@result, sum(A@values, B@values), 0.001)
omxCheckCloseEnough(model[['test17']]@result, prod(A@values, B@values), 0.001)
omxCheckCloseEnough(model[['test18']]@result, max(A@values, B@values), 0.001)
omxCheckCloseEnough(model[['test19']]@result, min(A@values, B@values), 0.001)
omxCheckCloseEnough(model[['test20']]@result, abs(A@values), 0.001)
omxCheckCloseEnough(model[['test21']]@result, cos(A@values), 0.001)
omxCheckCloseEnough(model[['test22']]@result, cosh(A@values), 0.001)
omxCheckCloseEnough(model[['test23']]@result, sin(A@values), 0.001)
omxCheckCloseEnough(model[['test24']]@result, sinh(A@values), 0.001)
omxCheckCloseEnough(model[['test25']]@result, tan(A@values), 0.001)
omxCheckCloseEnough(model[['test26']]@result, tanh(A@values), 0.001)
omxCheckCloseEnough(model[['test27']]@result, exp(A@values), 0.001)
omxCheckCloseEnough(model[['test28']]@result, log(A@values), 0.001)
omxCheckCloseEnough(model[['test29']]@result, sqrt(A@values), 0.001)
