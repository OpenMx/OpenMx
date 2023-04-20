#
#   Copyright 2007-2021 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


#------------------------------------------------------------------------------
# Author: Luis Araujo
# Date: 2022.09.28
# Filename: test-plus-overloading.R
# Purpose: Test if the plus operator works
#------------------------------------------------------------------------------


library(OpenMx)
library(testthat)
context("Plus operator")


mb <- mxModel("bivariate Heterogeneity Path Specification",
  type = "RAM",
  manifestVars = c("X", "Y")
) +
  mxPath(from = c("X", "Y"), arrows = 2, free = T, values = 1, lbound = .01) +
  mxPath(from = "X", to = "Y", arrows = 2, free = T, values = .2, lbound = .01) +
  mxPath(
    from = "one", to = c("X", "Y"), arrows = 1, free = T, values = c(0.1, -0.1),
    ubound = c(NA, 0), lbound = c(0, NA)
  )

mb <- mxGenerateData(mb, nrows = 1000, returnModel = T)
out <- try(mxRun(mb))

omxCheckEquals(is(out, "try-error"), FALSE)
