library(OpenMx)

omxCheckEquals(imxDetermineDefaultOptimizer(), "SLSQP")

omxCheckEquals(options()[['mxCondenseMatrixSlots']], FALSE)
