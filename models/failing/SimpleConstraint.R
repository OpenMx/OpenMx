# valgrind detects problems when run with CSOLNP

library(OpenMx)

#Create a constraint between MxMatrices 'A' and 'B'
constraint <- mxConstraint(A > B, name = 'AdominatesB')

# Constrain matrix 'K' to be equal to matrix 'limit'

model <- mxModel(model="con_test",
                 mxMatrix(type="Full", nrow=2, ncol=2, free=TRUE, name="K"),
                 mxMatrix(type="Full", nrow=2, ncol=2, free=FALSE, name="limit", values=
                            1:4),
                 mxConstraint(K == limit, name = "Klimit_equality"),
                 mxAlgebra(min(K), name="minK"),
                 mxAlgebraObjective("minK")
)
## Not run:

fit <- mxRun(model)
fit$matrices$K$values
