library(OpenMx)

# Rgui on Windows does not support stdout/stderr
#
# http://cran.r-project.org/bin/windows/base/rw-FAQ.html
# See "8.6 The output from my C code disappears. Why?"

# This should not result in an error
omxCheckTrue(imxLog("test"))
