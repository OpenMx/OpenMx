# stuff to dye for
# http://xxm.times.uh.edu/learn-xxm/lme4-example-dyestuff/

libraries <- rownames(installed.packages())
if (!("lme4" %in% libraries)) stop("SKIP")

library(lme4)
library(OpenMx)

batch <- mxModel(
  'batch', type="RAM",
  latentVars = c('batch'),
  mxPath('batch', arrows=2))

runTest <- function(data) {
  yield <- mxModel(
    'yield', type='RAM',
    mxModel(batch, data),
    manifestVars = c('Yield'),
    mxData(Dyestuff, 'raw'),
    mxPath('one', 'Yield'),
    mxPath('Yield', arrows=2, values=1, ubound=10000),
    mxPath('batch.batch', 'Yield', free=FALSE, values=1, joinKey="Batch"))
  yield <- mxRun(yield)
}

omxCheckError(runTest(mxData(data.frame(batch=unique(unclass(Dyestuff$Batch))),
         'raw')),
         "Attempt to join foreign key 'Batch' in yield.data of type 'unordered factor' with batch.data which has no primary key declared")

omxCheckError(runTest(mxData(data.frame(batch=unique(unclass(Dyestuff$Batch))),
         'raw', primaryKey='batch')),
         "Primary key 'batch' in batch.data of type 'integer' cannot be joined with foreign key 'Batch' in yield.data of type 'unordered factor'")

bad1 <- unique(Dyestuff$Batch)
levels(bad1)[3] <- 'Z'

omxCheckError(runTest(mxData(data.frame(batch=bad1),
                             'raw', primaryKey='batch')),
              "Primary key 'batch' in batch.data has different factor levels ('Z' != 'C') compared to foreign key 'Batch' in yield.data")

bad2 <- unique(Dyestuff$Batch)
levels(bad2) <- c(levels(bad2), 'Z')

omxCheckError(runTest(mxData(data.frame(batch=bad2),
                             'raw', primaryKey='batch')),
              "Primary key 'batch' in batch.data has a different number of factor levels compared to foreign key 'Batch' in yield.data")
