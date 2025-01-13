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

model1 <- mxModel(
	type = "RAM",
	manifestVars = c("x1", "x2"),
	mxPath(from = "x1", to = "x2", 
				 arrows = 1, 
				 labels = "b_x2_x1"),
	mxPath(from = "x1", 
				 to = "x2", 
				 arrows = 2, 
				 labels = "v_x2_x1",strictUnigraph=TRUE))
# adding the covariance removes the directed effect x1 -> x2
omxCheckEquals(names(coef(model1)),"v_x2_x1")
omxCheckTrue(all(coef(model1)==0))

model2 <- mxModel(
	type = "RAM",
	manifestVars = c("x1", "x2"),
	mxPath(from = "x1", to = "x2", 
				 arrows = 1, 
				 labels = "b_x2_x1"),
	mxPath(from = "x2", 
				 to = "x1", 
				 arrows = 2, 
				 labels = "v_x2_x1",strictUnigraph=TRUE))
# adding the covariance again remove the directed effect x1 -> x2
omxCheckEquals(names(coef(model2)),"v_x2_x1")
omxCheckTrue(all(coef(model2)==0))

# Now, use `strictUnigraph=FALSE` ####

model <- mxModel(
	type = "RAM",
	manifestVars = c("x1", "x2"),
	mxPath(
		from = "x1", to = "x2", 
		arrows = 1, 
		labels = "b_x2_x1"),
	mxPath(
		from = "x1", 
		to = "x2", 
		arrows = 2, 
		labels = "v_x2_x1",strictUnigraph=FALSE)
)
# Adding the covariance should not clobber the directed effect x1 -> x2
omxCheckEquals(names(coef(model)),c("b_x2_x1","v_x2_x1"))
omxCheckTrue(all(coef(model)==0))

model <- mxModel(
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
		labels = "b_x2_x1",strictUnigraph=FALSE)
)
# adding the directed effect should not clobber the covariance x1 <-> x2
omxCheckEquals(names(coef(model)),c("b_x2_x1","v_x2_x1"))
omxCheckTrue(all(coef(model)==0))
