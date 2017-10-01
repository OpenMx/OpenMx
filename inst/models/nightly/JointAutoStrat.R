libraries <- rownames(installed.packages())
if (!all(c("emx","gmp") %in% libraries)) stop("SKIP")

library(OpenMx)
library(emx)
library(rpf)
library(gmp)

numPatInclude <- c()
for (np in seq(2,200,2)) {
  numOutcome <- factorize(np)
  if (any(numOutcome > 7)) next
  numPatInclude <- c(numPatInclude, np)
}

todo <- expand.grid(numPat=numPatInclude, numContinuous=c(1L,100L),
                    numUnique=NA, strat=NA)

for (trial in 1:nrow(todo)) {
#for (trial in which(is.na(todo$continuous))) {
  set.seed(trial)
  numPat <- todo[trial,'numPat']
  exampleData <- list()
  numOutcome <- factorize(numPat)
  numOrdinal <- length(numOutcome)
  for (ix in 1:numOrdinal) {
    exampleData[[paste0('o',ix)]] <-
      mxFactor(1, levels=1:as.integer(numOutcome[ix]))
  }
  exampleData <- as.data.frame(exampleData)
  
  numContinuous <- todo[trial,'numContinuous']
  
  manifests <- c(paste0("o", 1:numOrdinal), paste0("c", 1:numContinuous))
  
  j1 <- mxModel(
    'j1', type="RAM",
    manifestVars=manifests,
    latentVars = 'G',
    emxThresholds(exampleData),
    mxData(exampleData, 'raw'),
    mxPath('G', manifests, values=runif(length(manifests), .5,1.5)),
    mxPath('G', free=FALSE, arrows=2, values=1),
    mxPath(paste0("o", 1:numOrdinal), arrows=2,
           values=.1+rlnorm(numOrdinal), free=FALSE),
    mxPath(paste0("c", 1:numContinuous), arrows=2,
           values=.1+rlnorm(numContinuous), free=TRUE),
    mxPath('one', paste0("c", 1:numContinuous),
           values=runif(numContinuous,-.5,.5)),
    mxComputeOnce('fitfunction','fit'))
  
  j1$expectation$thresholds <- "thresholdMatrix"
  
  ss <- 1000
  j2 <- mxGenerateData(j1, nrows=ss, returnModel=TRUE)
  cdf <- compressDataFrame(j2$data$observed[,paste0('o', 1:numOrdinal),drop=FALSE])
  todo[trial,'numUnique'] <- nrow(cdf)

  if (j2$fitfunction$jointConditionOn != 'auto') stop("not auto?")

#  j2$fitfunction$verbose <- 2L
  got <- mxRun(j2, silent = TRUE)

  strat <- attr(got$fitfunction$result,'jointConditionOn')
  if (!is.factor(todo$strat)) {
    todo$strat <- factor(todo$strat, levels=levels(strat))
  }
  todo[trial,c('strat')] <- strat
}

omxCheckTrue(todo$numUnique <= todo$numPat)
pred <- -3.58 + 0.36 * 1000/todo$numUnique - 0.06 * todo$numContinuous
omxCheckEquals(todo[pred<0, 'strat'], 'continuous')
omxCheckEquals(todo[pred>0, 'strat'], 'ordinal')
