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


# -----------------------------------------------------------------------
# Program: TestRowObjectiveExistenceVector.R  
#  Author: tbrick
#    Date: 2011-05-04
#
# Test of the OpenMx Row Objective Existence Vector
#
# Revision History
#     Timothy R. Brick -- Created
# -----------------------------------------------------------------------

require(OpenMx)
set.seed(200)
pctMissing=50
n = 200

# Generate data
# -----------------------------------------------------------------------
xes = cbind( rnorm(n, 0, 1), runif(n, 0, 1))
dimnames(xes) = list(NULL, c("x1", "x2"))
tSelect <- xes[,2] < (pctMissing/100)
xes[tSelect,2] <- NA
data = data.frame(xes)

# The objective is the count of non-missing values.
# -----------------------------------------------------------------------
sumModel = mxModel("sumModel",
	mxAlgebra(existenceVector, name="rowAlgebra"),
	mxAlgebra(sum(rowResults), name="reduceAlgebra"),
	mxData(data, type="raw"),
	mxFitFunctionRow(
		rowAlgebra='rowAlgebra',
		reduceAlgebra='reduceAlgebra',
		dimnames=c('x1','x2'),
	)
)
sumFit = mxRun(sumModel)
# The objective is the count of non-missing values.
# -----------------------------------------------------------------------
countModel = mxModel("countModel",
    mxMatrix("Full", 1, 2, free=FALSE, values=1, name="Unity"),
	mxAlgebra(sum(omxSelectCols(Unity, existenceVector)), name="rowAlgebra"),
	mxAlgebra(sum(rowResults), name="reduceAlgebra"),
	mxData(data, type="raw"),
	mxFitFunctionRow(
		rowAlgebra='rowAlgebra',
		reduceAlgebra='reduceAlgebra',
		dimnames=c('x1','x2'),
	)
)
countFit = mxRun(countModel)
# Correct answer: the count of all nonmissing values.
# -----------------------------------------------------------------------
realCount = sum(sum(!is.na(xes)))

omxCheckCloseEnough(realCount, as.vector(mxEval(objective, sumFit)))
omxCheckCloseEnough(!is.na(xes), mxEval(rowResults, sumFit))
omxCheckCloseEnough(realCount, as.vector(mxEval(objective, countFit)))
omxCheckCloseEnough(apply(!is.na(xes), 1, sum), as.vector(mxEval(rowResults, countFit)))
