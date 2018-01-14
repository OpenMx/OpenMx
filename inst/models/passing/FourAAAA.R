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

A <- mxMatrix(nrow = 1, ncol = 1, labels = 'data.A', name = 'A')
model1 <- mxModel('model1', A)

C <- mxAlgebra(model1.A, name = 'C')
model3 <- mxModel('model3', C)

B <- mxMatrix(nrow = 2, ncol = 1, free = TRUE, labels = c('A', 'data.A'), name = 'B')
data <- mxData(matrix(1, dimnames = list(c(), c('A'))), type = 'raw')
model2 <- mxModel('model2', B, model3, data)

superModel <- mxModel('superModel')
data <- mxData(matrix(0, dimnames = list(c(), c('A'))), type = 'raw')
superModel <- mxModel(superModel, model1, model2, data)

namespace <- imxGenerateNamespace(superModel)
flatModel <- imxFlattenModel(superModel, namespace)

omxCheckSetEquals(namespace$entities$superModel, c('model1', 'model2', 'data'))
omxCheckSetEquals(namespace$entities$model1, c('A'))
omxCheckSetEquals(namespace$entities$model2, c('B', 'data', 'model3'))
omxCheckSetEquals(namespace$entities$model3, c('C'))
