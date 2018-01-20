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
A <- mxMatrix('Full', 1, 1, name = 'A')
B <- mxAlgebra(A + A, name = 'B')
omxCheckError(mxModel('model', 'A', 'B'),
	paste("I don't know what to do with the following",
		"strings 'A' and 'B' that have been passed into the",
		"function: mxModel(\"model\", \"A\", \"B\")"))
model <- mxModel('model', A, B)
omxCheckError(mxModel(model, A, B, remove=TRUE),
	paste("Cannot use named entities",
	"when remove = TRUE. Instead give",
	"the name of the entity when removing it.",
	"See http://openmx.ssri.psu.edu/wiki/mxmodel-help#Remove_an_object_from_a_model"))
model <- mxModel(model, 'A', 'B', remove=TRUE)
omxCheckEquals(length(names(model)), 12)

genModel <- function() {
  mxModel("remove",
          mxMatrix(nrow=1, ncol=1, labels="x", name="X"),
          mxConstraint(X == 0, "e1"))
}

m1 <- genModel()
omxCheckTrue(!is.null(m1$e1))
m1$e1 <- NULL
omxCheckTrue(is.null(m1$e1))

m1 <- genModel()
m1[['e1']] <- NULL
omxCheckTrue(is.null(m1$e1))

m1 <- genModel()
m1 <- mxModel(m1, 'e1', remove=TRUE)
omxCheckTrue(is.null(m1$e1))
