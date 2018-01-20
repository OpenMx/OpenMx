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

omxCheckEquals(mxMakeNames(""), "i")
omxCheckEquals(mxMakeNames("103"), "X103")
omxCheckEquals(mxMakeNames("data"), "reserved")
omxCheckEquals(mxMakeNames("foo.bar[2,3]"), "fooxbarx2,3x")
omxCheckEquals(mxMakeNames("+!"), "xx")

omxCheckError(mxModel('abc -- def'), 
	paste("The name 'abc -- def' is illegal because",
	"it contains the '-' character in mxModel(\"abc -- def\")"))
omxCheckError(mxAlgebra(A + B, name = 'abc + 123|'),
	paste("The name 'abc + 123|' is illegal because",
	"it contains the characters '+' and '|' in",
	"mxAlgebra(A + B, name = \"abc + 123|\")"))
omxCheckError(mxMatrix('Full', 1, 1, labels = 'abc&def',
	name = 'B'), paste("The reference 'abc&def' in",
	"mxMatrix(\"Full\", 1, 1, labels = \"abc&def\", name = \"B\")",
	"is illegal because it contains the '&' character."))
omxCheckError(mxMatrix('Full', 1, 1, name = c('foo', 'bar', 'baz')),
	paste("The 'name' argument must be a single character",
	"argument in mxMatrix(\"Full\", 1, 1, name = c(\"foo\", \"bar\", \"baz\"))"))
omxCheckError(mxAlgebra(1 + 2 + 3, name = c('foo', 'bar', 'baz')),
	paste("The 'name' argument must be a single character",
	"argument in mxAlgebra(1 + 2 + 3, name = c(\"foo\", \"bar\", \"baz\"))"))
omxCheckError(mxConstraint(a == b, name = c('foo', 'bar', 'baz')),
	paste("The 'name' argument must be a single character",
	"argument in mxConstraint(a == b, name = c(\"foo\", \"bar\", \"baz\"))"))
omxCheckError(mxModel(name = c('foo', 'bar', 'baz')),
	paste("The 'name' argument must be a single character",
	"argument in mxModel(name = c(\"foo\", \"bar\", \"baz\"))"))
omxCheckError(mxMatrix('Full', 1, 1, labels = "foo.bar"),
	paste("The reference 'foo.bar' is illegal because it",
	"contains the '.' character in mxMatrix(\"Full\",",
	"1, 1, labels = \"foo.bar\") . To write a definition",
	"variable use 'data.bar'"))
omxCheckError(mxMatrix('Full', 1, 1, labels = "foo.bar.baz"),
	paste("The reference 'foo.bar.baz' is illegal because it",
	"contains the '.' character but it is not a valid definition",
	"variable in mxMatrix(\"Full\", 1, 1, labels = \"foo.bar.baz\")"))
omxCheckError(mxMatrix('Full', 1, 1, labels = "a.b.c.d.e"),
	paste("The reference 'a.b.c.d.e' is illegal because it",
	"contains the '.' character but it is not a valid definition",
	"variable in mxMatrix(\"Full\", 1, 1, labels = \"a.b.c.d.e\")"))

nameless <- mxMatrix('Full', name='abc', 1, 1)
nameless$name <- ''
omxCheckError(mxModel(model="complain", nameless),
	paste("Entity 'FullMatrix' in model 'complain' needs a name"))
