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
omxGraphviz(fit)

# dot has a manifest with cov arrows, but no means triangle
# digraph univSat1 { 
#    node [style=filled, fontname="Arial", fontsize=16]
#    X [shape=box, fillcolor="#a9fab1", height=0.5, width=0.5];
#    X -> X[dir=both, headport=s, tailport=s];
# }

# summary(fit)
#    name matrix row col   Estimate  Std.Error
# 1 var X      S   1   1 0.99325690 0.03114318
# 2  <NA>      M   1   1 0.01635556 0.02262332
# ...
# Observed statistics:  1000 
# Estimated parameters:  2 
