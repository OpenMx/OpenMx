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
a <- mxMatrix(values = matrix(1, 1, 1), name = 'a')
b <- mxMatrix(values = matrix(2, 1, 1), name = 'b')
c <- mxMatrix(values = matrix(3, 1, 1), name = 'c')
d <- mxAlgebra(a + b + c, name = 'd')

e <- mxMatrix(labels = c('a[1,1]','b[1,1]','c[1,1]','d[1,1]'), 
	nrow = 2, ncol = 2, byrow = TRUE, name = 'e')

model <- mxModel(a, b, c, d, e)

model <- mxRun(model)

expected <- matrix(c(1, 2, 3, 6), nrow = 2, ncol = 2, byrow = TRUE)

omxCheckCloseEnough(model[['e']]@values, expected, 0.0001)
