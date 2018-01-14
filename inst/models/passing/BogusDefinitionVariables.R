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

model <- mxModel('model')
foo <- mxMatrix('Full', 1, 1, labels = 'bar.data.baz', name = 'foo')
model <- mxModel(model, foo)
omxCheckError(mxRun(model), paste("The definition variable",
	"'bar.data.baz' in matrix 'foo'",
	"refers to a data set that does not exist"))
foo <- mxMatrix('Full', 1, 1, labels = 'data.baz', name = 'foo')
model <- mxModel(model, foo)
omxCheckError(mxRun(model), paste("A definition variable 'data.baz'",
	"has been declared in model 'model'",
	"that does not contain a data set"))
data <- mxData(matrix(1,1,1), type = 'raw')
model <- mxModel(model, data)
omxCheckError(mxRun(model), paste("The definition variable",
	"'model.data.baz' in matrix 'foo'",
	"refers to a data set that does not contain column names"))
data <- mxData(matrix(1,1,1, dimnames = list('a', 'b')), type = 'raw')
model <- mxModel(model, data)
omxCheckError(mxRun(model), paste("The definition variable",
	"'model.data.baz' in matrix 'foo'",
	"refers to a data set that does",
	"not contain a column with name 'baz'"))
