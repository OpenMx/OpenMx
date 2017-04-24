library(OpenMx)

test1 <- mxModel("test1",
	mxMatrix(values=1, nrow=1, ncol=1, free=TRUE, name="thang"),
	mxMatrix(nrow=1, ncol=1, labels="abscissa1", free=TRUE, name="currentAbscissa"),
	mxMatrix(values=-2:2, nrow=1, ncol=5, name="abscissa",
		 dimnames=list(c('abscissa1'), NULL)),
	mxAlgebra(cbind(currentAbscissa + thang, currentAbscissa * thang), name="stuff"),
	mxAlgebra(mxEvaluateOnGrid(stuff, abscissa), name="grid"))
	
omxCheckError(mxRun(test1), "mxEvaluateOnGrid: algebra 'test1.stuff' returned 2 columns instead of 1");

test2 <- mxModel("test2",
	mxMatrix(values=1.1, nrow=1, ncol=1, free=TRUE, name="thang"),
	mxMatrix(nrow=1, ncol=1, labels="abscissa1", free=TRUE, name="currentAbscissa"),
	mxMatrix(values=-2:2, nrow=1, ncol=5, name="abscissa",
		 dimnames=list(c('abscissa1'), NULL)),
	mxAlgebra(rbind(currentAbscissa + thang, currentAbscissa * thang), name="stuff"),
	mxAlgebra(mxEvaluateOnGrid(stuff, abscissa), name="grid"))
	
test2 <- mxRun(test2)
omxCheckCloseEnough(test2$grid$result, matrix(c(-1:3 + .1, -2:2 * 1.1), ncol=5, nrow=2,byrow=TRUE))

test3 <- mxModel("test3",
	mxMatrix(values=1.1, nrow=1, ncol=1, free=TRUE, name="thang"),
	mxMatrix(nrow=1, ncol=1, labels="abscissa1", free=TRUE, name="currentAbscissa"),
	mxMatrix(values=-2:2, nrow=1, ncol=5, name="abscissa"),
	mxAlgebra(rbind(currentAbscissa + thang, currentAbscissa * thang), name="stuff"),
	mxAlgebra(mxEvaluateOnGrid(stuff, abscissa), name="grid"))

omxCheckError(mxRun(test3), "mxEvaluateOnGrid: abscissa 'test3.abscissa' must have rownames")

test4 <- mxModel("test4",
	mxMatrix(values=1.1, nrow=1, ncol=1, free=TRUE, name="thang"),
	mxMatrix(nrow=1, ncol=1, labels="abscissa1", free=TRUE, name="currentAbscissa"),
	mxMatrix(values=-2:2, nrow=1, ncol=5, name="abscissa",
		 dimnames=list(c('xyz'), NULL)),
	mxAlgebra(rbind(currentAbscissa + thang, currentAbscissa * thang), name="stuff"),
	mxAlgebra(mxEvaluateOnGrid(stuff, abscissa), name="grid"))

omxCheckError(mxRun(test4), "mxEvaluateOnGrid: abscissa 'test4.abscissa' row 1, 'xyz' does not name a free parameter")
