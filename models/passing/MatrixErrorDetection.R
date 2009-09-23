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
A <- mxMatrix('Full', 1, 1, labels = 'data.foo', free = TRUE, name = 'A')
model <- mxModel('model', A)
omxCheckError(mxRun(model), 
	paste("The definition variable 'data.foo'",
		"has been assigned to a",
		"free parameter in matrix 'A'"))
A <- mxMatrix('Full', 1, 1, labels = 'foo[1,2]', free = TRUE, name = 'A')
model <- mxModel('model', A)
omxCheckError(mxRun(model), 
	paste("The substitution 'foo[1,2]'",
		"has been assigned to a",
		"free parameter in matrix 'A'"))
