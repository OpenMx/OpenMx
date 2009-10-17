# SCRIPT: omx_graphvis_example.R
# Author: Timothy  Bates tim.bates@ed.ac.uk
# History:  Fri Oct 16 17:02:32 BST 2009
# OpenMx: http://www.openmx.virginia.com
#  .dot output not modelling means
##########################################
require(OpenMx)
testData = as.matrix(rnorm(n=1000, mean=0, sd=1)) # make a thousand numbers, with mean 0 and sd 1
selVars = c("X")
colnames(testData) <- selVars
model <- mxModel("univSat1", manifestVars= selVars,
  mxPath(from="X", arrows=2, free=TRUE, values=1, lbound=.01, labels="var X"), 
  mxData(observed=testData, type="raw"), 
  type="RAM"
)
fit = mxRun(model)
omxGraphviz(fit,'test.dot')
# system("open test.dot") # uncomment this to have R open the dot file in your default app (omnigraffle etc)
# getwd() # will show you where the file is written