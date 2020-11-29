#
#   Copyright 2007-2018 by the individuals mentioned in the source code history
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

library(OpenMx)
library(testthat)

if(mxOption(key="Default optimizer") != 'SLSQP') stop("SKIP")

varNames <- c('x','y','z')

origCov <- matrix(c(3.6,2.2,0.5,2.2,3.2,2,0.5,2,3.9), nrow = 3,
	dimnames = list(varNames, varNames))

data <- mxData(origCov, type = "cov", numObs = 10)

s <- mxMatrix("Symm", free = TRUE,
	values = matrix(c(10,9,9,9,10,9,9,9,10), nrow = 3),
	labels = matrix(c("v1", "c12", "c13",
		"c12", "v2", "c23",
		"c13", "c23", "v3"), nrow = 3),
	dimnames = list(varNames, varNames),
	name = "s", lbound = -100, ubound = 100)

c <- mxMatrix("Full", free = FALSE, values = 0.4, nrow = 3, ncol=3, name= "c")

model <- mxModel("model", data, s, c,
		 mxConstraint(s > c, name = "constraint"),
		 mxBounds(c("v1", "v2", "v3"), min = 0),
		 mxFitFunctionML(),
		 mxExpectationNormal("s"))

run <- mxRun(model)

omxCheckCloseEnough((mxEval(s, run)/ origCov), matrix(.9,3,3), 1e-4) # constraint not active

c <- mxMatrix("Full", free = FALSE, values = 0.6, nrow = 3, ncol=3, name= "c")

model <- mxModel("model", data, s, c,
		 mxConstraint(c < s, name = "constraint"),
		 mxBounds(c("v1", "v2", "v3"), min = 0),
		 mxFitFunctionML(),
		 mxExpectationNormal("s"))

# Constraint is duplicated so tolerance is halved
#model <- mxOption(model, key="feasibility tolerance", value = .05)

# No code 6 warning because constraints are active
run <- expect_warning(mxRun(model), NA)

# gradient should not be too close to zero
expect_equal(max(abs(run$output$gradient)), .267, .01)

omxCheckTrue(coef(run)['c13'] < 0.6)  # constraint should be active!
omxCheckCloseEnough(var(c(mxEval(s, run) / origCov)), 0, 0.011)
#cat(deparse(round(c(vech(mxEval(s, run))),3)))
omxCheckCloseEnough(vech(mxEval(s, run)), c(3.284, 2.066, 0.581, 2.979, 1.892, 3.56), .001)
