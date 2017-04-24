# ===========
# = History =
# ===========
# 2017-04-14 04:50PM TBATES Correct bug in error check - this is now passing.

#
#   Copyright 2007-2012 The OpenMx Project
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
# 2013-01-10 update for new error reporting when using named entity to constrain label


require(OpenMx)
# Check we error when a free named entity shares its name with a label in an object
foo <- mxMatrix(name = 'foo', nrow = 1, ncol = 1, free = TRUE, labels = 'foo')
bar <- mxAlgebra(foo, name = 'bar')
model <- mxModel('model', foo, bar)
expErr = "In model 'model' the following are both named entities and free parameters: 'foo' 
If you are trying to set a path using an mxAlgebra, then refer to the Algebra with square-bracket notation. 
i.,e, instead of labels=\" 'foo' \" use: labels=\" 'foo' [1,1]\""

omxCheckError(mxRun(model), expErr)

# paste("In model 'model' the following are both named entities and free parameters: 'foo'",
# "\nIf you are trying to set a path using an mxAlgebra, then refer to the Algebra with square-bracket notation.",
# "\ni.,e, instead of labels=\"", omxQuotes(overlap), "\" use: labels=\"", omxQuotes(overlap), "[1,1]\"")

foo   <- mxMatrix(name = 'foo',  nrow = 1, ncol = 1, free = TRUE , labels = 'a')
bar   <- mxMatrix(name = 'bar',  nrow = 1, ncol = 1, free = TRUE , labels = 'a')
baz   <- mxMatrix(name = 'baz',  nrow = 1, ncol = 1, free = FALSE, labels = 'a')
quux  <- mxMatrix(name = 'quux', nrow = 1, ncol = 1, free = FALSE, labels = 'a')
model <- mxModel('model', foo, bar, baz, quux)
omxCheckError(mxRun(model), "In model 'model' the name 'a' is used as a free parameter in 'model.foo' and 'model.bar' and as a fixed parameter in 'model.baz' and 'model.quux'")
