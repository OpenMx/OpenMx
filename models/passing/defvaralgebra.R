require(OpenMx)

A <- mxMatrix('Full', 1, 1, labels = 'data.foo', name = 'A')
B <- mxAlgebra(data.foo + A, name = 'B')
data <- mxData(matrix(c(1:5), 5, 1, dimnames = list(NULL, c('foo'))), type='raw')
model <- mxModel('model', A, B, data)

omxCheckEquals(mxEval(A[1,1], model, compute=FALSE), 0)
omxCheckError(mxEval(A[1,1], model, compute=TRUE, defvar.row = 0),
	paste("Row number '0' is out of",
	"bounds for definition variable",
	"'data.foo' used in the context of 'A'\nOne possibility is you didn't add this variable to the data for this group?"))
omxCheckEquals(mxEval(A[1,1], model, compute=TRUE), 1)
omxCheckEquals(mxEval(A[1,1], model, compute=TRUE, defvar.row=2), 2)
omxCheckEquals(mxEval(A[1,1], model, compute=TRUE, defvar.row=3), 3)
omxCheckEquals(mxEval(A[1,1], model, compute=TRUE, defvar.row=4), 4)
omxCheckEquals(mxEval(A[1,1], model, compute=TRUE, defvar.row=5), 5)
omxCheckError(mxEval(A[1,1], model, compute=TRUE, defvar.row = 6),
	paste("Row number '6' is out of",
	"bounds for definition variable",
	"'data.foo' used in the context of 'A'\nOne possibility is you didn't add this variable to the data for this group?"))

omxCheckEquals(mxEval(B[1,1], model, compute=TRUE), 2)
omxCheckEquals(mxEval(B[1,1], model, compute=TRUE, defvar.row=2), 4)
omxCheckEquals(mxEval(B[1,1], model, compute=TRUE, defvar.row=3), 6)
omxCheckEquals(mxEval(B[1,1], model, compute=TRUE, defvar.row=4), 8)
omxCheckEquals(mxEval(B[1,1], model, compute=TRUE, defvar.row=5), 10)
omxCheckError(mxEval(B[1,1], model, compute=TRUE, defvar.row = 0),
	paste("Row number '0' is out of",
	"bounds for definition variable",
	"'data.foo' used in the context of 'B'\nOne possibility is you didn't add this variable to the data for this group?"))
omxCheckError(mxEval(B[1,1], model, compute=TRUE, defvar.row = 6),
	paste("Row number '6' is out of",
	"bounds for definition variable",
	"'data.foo' used in the context of 'B'\nOne possibility is you didn't add this variable to the data for this group?"))
