library(OpenMx)

data(multiData1)
manifests <- c("x1", "y")
multiData1Cov <- cov(multiData1[,c(1,5)])

uniRegModel <- mxModel("Univariate Regression of y on x1",
    type="RAM",
    manifestVars=manifests,
    mxPath(from="x1", to="y", arrows=1, 
           free=TRUE, values=.2, labels="b1"),
    mxPath(from=manifests, arrows=2, 
           free=TRUE, values=-.8, labels=c("VarX1", "VarE")),
    mxData(observed=multiData1Cov, type="cov", numObs=500)
    )

ign <- omxCheckWarning(try(mxRun(uniRegModel)),
                       paste("In model 'Univariate Regression of y on x1' Optimizer returned a non-zero status code 10.",
                             "Starting values are not feasible. Consider mxTryHard()"))

# ---------- Newton Raphson

mat1 <- mxMatrix("Full", rnorm(1), free=TRUE, nrow=1, ncol=1, labels="m1", name="mat1")
obj <- mxAlgebra(1/0, name = "obj")
grad <- mxAlgebra(1, name = "grad", dimnames=list("m1", NULL))
hess <- mxAlgebra(.5, name = "hess", dimnames=list("m1", "m1"))

model1 <- mxModel("model1", mat1, obj, grad, hess,
                  mxFitFunctionAlgebra("obj", gradient="grad", hessian="hess"),
                  mxComputeNewtonRaphson())
ign <- omxCheckWarning(mxRun(model1), paste("In model 'model1' Optimizer returned a non-zero status code 10.",
                                            "Starting values are not feasible. Consider mxTryHard()"))
