# check to ensure that the default optimizer is set to NPSOL
omxCheckEquals(imxDetermineDefaultOptimizer(), "NPSOL")

omxCheckEquals(options()[['mxCondenseMatrixSlots']], FALSE)
