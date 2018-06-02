library(OpenMx)

m1 <- mxModel(
	"colville",
	mxMatrix(type='Full', values=runif(4, min=-10,max=10),
		nrow=4, ncol=1, free=TRUE, lbound=-10, ubound=10, name='x'),
	mxAlgebra(100.*(x[2,1]-(x[1,1]*x[1,1]))*(x[2,1]-(x[1,1]*x[1,1]))
		+ (1.-x[1,1])*(1.-x[1,1])
		+ 90.*(x[4,1]-(x[3,1]*x[3,1]))*(x[4,1]-(x[3,1]*x[3,1]))
		+ (1.-x[3,1])*(1.-x[3,1])
		+ 10.1*((x[2,1]-1.)*(x[2,1]-1.) + (x[4,1]-1.)*(x[4,1]-1.))
		+ 19.8*(x[2,1]-1.)*(x[4,1]-1.), name="fit"),
	mxFitFunctionAlgebra('fit'),
	mxComputeSimAnnealing(method='ingber2012',
		control=list('Limit_Acceptances'=50000,
			'Temperature_Anneal_Scale'=1000,
			'Cost_Parameter_Scale_Ratio'=1.0E-1,
			'Accepted_To_Generated_Ratio'=1.0E-4,
			'Limit_Invalid_Generated_States'=1e9
		)))
m1 <- mxRun(m1)
omxCheckCloseEnough(m1$output$fit, 0, 1e-2)
omxCheckCloseEnough(coef(m1), rep(1,4), .1)

# ROSENBROCK, H. H. An automatic method for finding the greatest or
# least value of a function.  Comput. J. 3 (1960), 175-184.

m1 <- mxModel(
	"rosenbrock2",
	mxMatrix(type='Full', values=c(-1.2, 1),
		nrow=2, ncol=1, free=TRUE, lbound=-2.048, ubound=2.048, name='x'),
	mxAlgebra(100*(x[2,1]-x[1,1]^2)^2 + (1-x[1,1])^2, name="fit"),
	mxFitFunctionAlgebra('fit'),
	mxComputeSimAnnealing(verbose=0L, method='ingber2012',
		control=list(
			'Limit_Acceptances'=1e5,
			'Limit_Generated'=1e6,
			'Temperature_Ratio_Scale'=1e-2,
			'Initial_Parameter_Temperature'=1e10,
			'Maximum_Cost_Repeat'=0,
			'Reanneal_Cost'=5,
			'Reanneal_Parameters'=0
		)))

m1 <- mxRun(m1)
omxCheckCloseEnough(m1$output$fit, 0, 1e-2)
omxCheckCloseEnough(coef(m1), rep(1,2), .1)

set.seed(1234)
dimension <- 8

# Rastrgin function has global minimum at 0 with all parameters at 0.

m1 <- mxModel(
	"rastrigin",
	mxMatrix(type='Full', values=runif(dimension, min=-1e6, max=1e6),
		ncol=1, nrow=dimension, free=TRUE, name='x'),
	mxAlgebra(sum(x*x - 10 * cos(2 * pi * x)) + 10 * dimension, name="fit"),
	mxFitFunctionAlgebra('fit'),
	mxComputeSimAnnealing())

m1 <- mxRun(m1)

omxCheckCloseEnough(m1$output$fit, 0, 1e-2)
omxCheckCloseEnough(coef(m1), rep(0,dimension), 1e-2)

# ----

path <- paste0(tempdir(),'/','log.csv')

m2 <- mxModel(
	"rastrigin",
	mxMatrix(type='Full', values=runif(dimension, min=-1e6, max=1e6),
		ncol=1, nrow=dimension, free=TRUE, name='x'),
	mxAlgebra(sum(x*x - 10 * cos(2 * pi * x)) + 10 * dimension, name="fit"),
	mxFitFunctionAlgebra('fit'),
	mxComputeCheckpoint(path=path, append=FALSE))

m2 <- mxRun(m2)
m2 <- mxRun(m2)

log1 <- read.table(path, header=TRUE, sep='\t', check.names=FALSE)
omxCheckEquals(nrow(log1), 1)

m2$compute$append <- TRUE
m2$compute$header <- FALSE

m2 <- mxRun(m2)

log2 <- read.table(path, header=TRUE, sep='\t', check.names=FALSE)
omxCheckEquals(nrow(log2), 2)

# ----

dimension <- 2

m1 <- mxModel(
	"rastrigin",
	mxMatrix(type='Full', values=runif(dimension, min=-1e6, max=1e6),
		ncol=1, nrow=dimension, free=TRUE, name='x'),
	mxAlgebra(sum(x*x - 10 * cos(2 * pi * x)) + 10 * dimension, name="fit"),
	mxConstraint(x[1,1] > .5, "con1"),
	mxFitFunctionAlgebra('fit'),
	mxComputeSimAnnealing()
)

m1 <- mxRun(m1)
omxCheckCloseEnough(coef(m1), c(1,0), 1e-2)
