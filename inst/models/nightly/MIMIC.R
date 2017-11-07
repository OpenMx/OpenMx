# From http://openmx.ssri.psu.edu/wiki/example-models#MIMIC
require(OpenMx)

data = structure(c(1, 0.304, 0.305, 0.1, 0.284, 0.176, 0, 1, 0.344,
                   0.156, 0.192, 0.136, 0, 0, 1, 0.158, 0.324, 0.226, 0, 0, 0, 1, 
                   0.36, 0.21, 0, 0, 0, 0, 1, 0.265, 0, 0, 0, 0, 0, 1),
                 .Dim = c(6L,  6L),
                 .Dimnames = list(c("X1", "X2", "X3", "X4", "X5", "X6"), 
                                  c("X1", "X2", "X3", "X4", "X5", "X6")))
 
# now letsfill in the upper triangle with a flipped version of the lower
data[upper.tri(data, diag=F)] = t(data)[upper.tri(data, diag=F)]
 
# Set up manifest variables
manifests = c("income", "occup", "educ", "church", "member", "friends")
 
# Use these to create names for our dataframe
dimnames(data) = list(manifests, manifests)
 
# And latents
latents   = "social" # 1 latent, with three formative inputs, and three reflective outputs (each with residuals)
 
# Just to be helpful to myself, I've made lists of the formative sources, and the reflective receiver variables in this MIMIC model
receivers = manifests[4:6]
sources   = manifests[1:3]
 
MIMIC <- mxModel("MIMIC", type="RAM",
    manifestVars = manifests,
    latentVars   = latents,
 
    # Factor loadings
    mxPath(from = sources , to = "social", lbound=0),
    mxPath(from = "social", to = receivers, lbound=0),
 
    # Correlated formative sources for F1, each with variance = 1
    mxPath(from = sources, connect = "unique.bivariate", arrows = 2),
    mxPath(from = sources, arrows = 2, values = 1, free=F ),
 
    # Residual variance on receivers
    mxPath(from = receivers, arrows = 2),
    mxData(data, type = "cov", numObs = 530)
)

MIMIC <- mxRun(MIMIC)

omxCheckCloseEnough(MIMIC$output$fit, 2893.2189, .01)

got <- mxCheckIdentification(MIMIC)
omxCheckEquals(match(got$non_identified_parameters, names(MIMIC$output$estimate)), 1:6)
