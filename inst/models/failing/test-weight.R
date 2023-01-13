library(OpenMx)
library(testthat)

data("nhanesDemo")


m <- mxModel("testing weitghs", type = "RAM", manifestVars =c("fpl", "age", "psu"),
  mxPath(from = "age", to = "fpl", arrows = 1, free = T, values = 0.5),
  mxPath(from = "psu", to = "fpl", arrows = 1, free = T, values = 0.5),
  mxPath(from = "one", to = c("fpl", "age", "psu"), arrows = 1, free = T, values = 0.5),
  mxPath(from = "age", to = "age", arrows = 2, free = T, values = 1),
  mxPath(from = "psu", to = "psu", arrows = 2, free = T, values = 1),
  mxPath(from = "fpl", to = "fpl", arrows = 2, free = T, values = 1),
  mxData(nhanesDemo, type = "raw", weight = "persWeight")
 )

mo <- try(mxRun(m))

omxCheckEquals(is(mo, "try-error"), FALSE)
