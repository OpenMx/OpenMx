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
omxCheckError(mxModel('abc -- def'), 
	paste("The name 'abc -- def' is illegal because",
	"it contains the '-' character."))
omxCheckError(mxAlgebra(A + B, name = 'abc + 123|'),
	paste("The name 'abc + 123|' is illegal because",
	"it contains the characters '+' and '|'."))
omxCheckError(mxMatrix('Full', 1, 1, labels = 'abc&def',
	name = 'B'), paste("The reference 'abc&def' in matrix",
	"'B' is illegal because it contains the '&' character."))