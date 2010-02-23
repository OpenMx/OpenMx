#
#   Copyright 2007-2010 The OpenMx Project
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
data <- mxData(type = 'raw', matrix(".", 3, 3, dimnames = list(NULL,c('a','b','c'))))
covariance <- mxMatrix('Symm', 3, 3, values = c(1:6), name = 'cov')
means <- mxMatrix('Full', 1, 3, values = c(1:3), name = 'means')
objective <- mxFIMLObjective('cov', 'means')
model <- mxModel('model', objective, covariance, means, data)
omxCheckError(mxRun(model), paste("The data object", omxQuotes("model.data"),
	"contains an observed matrix that is not of type 'double'"))
