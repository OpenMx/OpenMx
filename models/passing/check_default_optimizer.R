# check to ensure that the default optimizer is set to NPSOL
omxCheckEquals(mxOption(NULL, "Default optimizer"), "NPSOL")
