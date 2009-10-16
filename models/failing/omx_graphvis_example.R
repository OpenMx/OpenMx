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
# Error: mxPath() call will generate 0 paths but you have specified 1 arrows with 'from' argument assigned to  and 'to' argument assigned to 
