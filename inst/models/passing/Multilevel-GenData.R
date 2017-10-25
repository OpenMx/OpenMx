library(OpenMx)

set.seed(1)

df <- NULL
for (batch in 1:50) {
  df <- rbind(df, expand.grid(case=1:(5+sample.int(4, 1)), Batch=batch, Yield=0))
}

batch <- mxModel(
    'batch', type="RAM",
    latentVars = c('batch'),
    mxData(data.frame(batch=unique(df$Batch)), 'raw', primaryKey='batch'),
    mxPath('batch', arrows=2, values=.75, lbound=.001))

trueYield <- mxModel(
    'yield', type='RAM', batch,
    manifestVars = c('Yield'),
    mxData(df, 'raw'),
    mxPath('one', 'Yield', values=1e-6),
    mxPath('Yield', arrows=2, values=1),
    mxPath('batch.batch', 'Yield', free=FALSE, values=1, joinKey="Batch"))

result <- expand.grid(rep=1:20)
for (px in names(coef(trueYield))) result[[px]] <- NA
result$rep <- NULL

for (rep in 1:nrow(result)) {
  yield <- mxGenerateData(trueYield, returnModel = TRUE)
  
  yield <- mxRun(yield, silent = TRUE)
  
  result[rep, names(coef(yield))] <- coef(yield)
}

omxCheckCloseEnough(colMeans(result) - coef(trueYield), rep(0,3), .06)
omxCheckCloseEnough(apply(result, 2, var), rep(0,3), .03)
