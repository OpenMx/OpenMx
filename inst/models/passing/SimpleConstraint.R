library(OpenMx)

#mxOption(NULL, "Default optimizer", "NPSOL")

# Constrain matrix 'K' to be equal to matrix 'limit'
model <- mxModel(model="con_test",
                 mxMatrix(type="Full", nrow=2, ncol=2, free=TRUE, name="K"),
                 mxMatrix(type="Full", nrow=2, ncol=2, free=FALSE, name="limit", values=1:4),
                 mxConstraint(K == limit, name = "Klimit_equality"),
                 mxAlgebra(min(K), name="minK"),
                 mxFitFunctionAlgebra("minK")
)

fit <- mxRun(model)

omxCheckCloseEnough(fit$output$fit, 1, 1e-4)
omxCheckCloseEnough(c(fit$matrices$K$values), 1:4, 1e-3)
