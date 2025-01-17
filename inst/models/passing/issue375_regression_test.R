#
#   Copyright 2007-2024 by the individuals mentioned in the source code history
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

# Adapted from syntax by Jannik Orzek & Michael D. Hunter

library(OpenMx)
# This script doesn't test anything optimization-related, so it doesn't need 
# to be run with anything other than the on-load default:
if(mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")


msg <- paste(
	"Looks like there is a pre-existing one-headed path from 'x1' to 'x2'.\n",
	"That path is now overwritten by a two-headed path between 'x1' and 'x2'.\n",
	"To retain the one-headed path, either use 'dummy' latent variables, or directly modify the MxModel's 'A' matrix;\n",
	"See the `mxPath()` help page for examples.\n",
	"Be advised, this overwriting behavior may change in the future, so do not write scripts that rely upon it!",
	sep=""
)
omxCheckWarning(
	model1 <- mxModel(
		type = "RAM",
		manifestVars = c("x1", "x2"),
		mxPath(from = "x1", to = "x2", 
					 arrows = 1, 
					 labels = "b_x2_x1"),
		mxPath(from = "x1", 
					 to = "x2", 
					 arrows = 2, 
					 labels = "v_x2_x1")),
	msg
)
# adding the covariance removes the directed effect x1 -> x2
omxCheckEquals(names(coef(model1)),"v_x2_x1")
omxCheckTrue(all(coef(model1)==0))


omxCheckWarning(
	model2 <- mxModel(
		type = "RAM",
		manifestVars = c("x1", "x2"),
		mxPath(from = "x1", to = "x2", 
					 arrows = 1, 
					 labels = "b_x2_x1"),
		# Note that x2 & x1 are swapped here (https://github.com/OpenMx/OpenMx/issues/375#issuecomment-2582989383):
		mxPath(from = "x2", 
					 to = "x1", 
					 arrows = 2, 
					 labels = "v_x2_x1")),
	msg
)
# adding the covariance removes the directed effect x1 -> x2
omxCheckEquals(names(coef(model2)),"v_x2_x1")
omxCheckTrue(all(coef(model2)==0))


msg <- paste(
	"Looks like there is a pre-existing two-headed path between 'x1' and 'x2'.\n",
	"That path is now overwritten by a one-headed path from 'x1' to 'x2'.\n",
	"To retain the two-headed path, either use 'dummy' latent variables, or directly modify the MxModel's 'S' matrix;\n",
	"See the `mxPath()` help page for examples.\n",
	"Be advised, this overwriting behavior may change in the future, so do not write scripts that rely upon it!",
	sep=""
)
omxCheckWarning(
	model3 <- mxModel(
	type = "RAM",
	manifestVars = c("x1", "x2"),
	mxPath(
		from = "x1", 
		to = "x2", 
		arrows = 2, 
		labels = "v_x2_x1"),
	mxPath(
		from = "x1", to = "x2", 
		arrows = 1, 
		labels = "b_x2_x1")),
	msg
)
omxCheckEquals(names(coef(model3)),"b_x2_x1")
omxCheckTrue(all(coef(model3)==0))


# Warning message should be the same as that for `model3`, due to symmetry of 'S' matrix:
omxCheckWarning(
	model4 <- mxModel(
		type = "RAM",
		manifestVars = c("x1", "x2"),
		# Again, note swap:
		mxPath(
			from = "x2", 
			to = "x1", 
			arrows = 2, 
			labels = "v_x2_x1"),
		mxPath(
			from = "x1", to = "x2", 
			arrows = 1, 
			labels = "b_x2_x1")),
	msg
)
omxCheckEquals(names(coef(model4)),"b_x2_x1")
omxCheckTrue(all(coef(model4)==0))


msg <- paste(
	"Looks like there is a pre-existing two-headed path between 'x2' and 'x1'.\n",
	"That path is now overwritten by a one-headed path from 'x2' to 'x1'.\n",
	"To retain the two-headed path, either use 'dummy' latent variables, or directly modify the MxModel's 'S' matrix;\n",
	"See the `mxPath()` help page for examples.\n",
	"Be advised, this overwriting behavior may change in the future, so do not write scripts that rely upon it!",
	sep=""
)
omxCheckWarning(
	model5 <- mxModel(
		type = "RAM",
		manifestVars = c("x1", "x2"),
		mxPath(
			from = "x1", 
			to = "x2", 
			arrows = 2, 
			labels = "v_x2_x1"),
		mxPath(
			from = "x2", to = "x1", 
			arrows = 1, 
			labels = "b_x1_x2")),
	msg
)
omxCheckEquals(names(coef(model5)),"b_x1_x2")
omxCheckTrue(all(coef(model4)==0))
