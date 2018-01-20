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
bar <- matrix(c(1:4), 2, 2)
foo <- mxAlgebra((-1 + 2.0 + 5e-1) %x% bar, 'foo')
model <- mxModel('model', foo)
modelOut <- mxRun(model)
omxCheckIdentical(matrix(c(1.5, 3.0, 4.5, 6.0), 2, 2), mxEval(foo, modelOut))
