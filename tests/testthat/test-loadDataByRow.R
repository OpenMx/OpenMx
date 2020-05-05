library(OpenMx)
library(testthat)
context("loadDataByRow")
suppressWarnings(RNGversion("3.5"))
set.seed(1)

if (mxOption(NULL,"Default optimizer") == 'NPSOL') stop("SKIP")
#mxOption(NULL, "Number of Threads", 1L)

data("jointdata", package ="OpenMx")

# specify ordinal columns as ordered factors
jointdata[,c(2,4,5)] <- mxFactor(jointdata[,c(2,4,5)], 
	levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

satCov <- mxMatrix("Symm", 5, 5,
	free=TRUE, values=diag(5), name="C")
satCov$free[2,2] <- FALSE
satCov$free[4,4] <- FALSE
satCov$free[5,5] <- FALSE

loadings <- mxMatrix("Full", 1, 5,
	free=TRUE, values=1, name="L", lbound=0)
loadings$ubound[1,4:5] <- 2

resid <- mxMatrix("Diag", 5, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=.5, name="U")
	
means <- mxMatrix("Full", 1, 5,
	free=c(TRUE, FALSE, TRUE, FALSE, FALSE), values=0, name="M")
	
thresh <- mxMatrix("Full", 3, 3, FALSE, 0, name="T")

thresh$free[,1] <- c(TRUE, FALSE, FALSE)
thresh$values[,1] <- c(0, NA, NA)
thresh$labels[,1] <- c("z2t1", NA, NA)

thresh$free[,2] <- TRUE
thresh$values[,2] <- c(-1, 0, 1)
thresh$labels[,2] <- c("z4t1", "z4t2", "z4t3")

thresh$free[,3] <- c(TRUE, TRUE, FALSE)
thresh$values[,3] <- c(-1, 1, NA)
thresh$labels[,3] <- c("z5t1", "z5t2", NA)

model1 <- mxModel("loadData",
	loadings, resid, means, thresh,
	mxAlgebra(t(L) %*% L + U, name="C"),
	mxFitFunctionWLS(),
	mxExpectationNormal("C", "M",
		dimnames=names(jointdata)[1:5],
		thresholds="T",
		threshnames=c("z2", "z4", "z5")))

result1 <- c()
numSets <- 8
dsets <- list()
for (dx in 1:numSets) {
  df <- data.frame(z1=jointdata$z1 + rnorm(nrow(jointdata), mean=dx/10, sd=.1),
                   z2=ordered(sample.int(2, nrow(jointdata), replace=TRUE)-1L))
  df$z1[sample.int(nrow(jointdata), sample(20,1))] <- NA
  df$z2[sample.int(nrow(jointdata), sample(20,1))] <- NA
  dsets[[dx]] <- df
  for (cx in paste0('z',3:5)) df[[cx]] <- jointdata[[cx]]
  model2 <- mxModel(model1, mxData(df, 'raw'), mxFitFunctionWLS())
  model2 <- mxRun(model2)
  result1 <- rbind(result1, c(coef(model2), model2$output$standardErrors))
}
flat <- t(sapply(unlist(dsets, recursive=FALSE), as.character))
colnames(flat) <- paste0("p",1:ncol(flat))
rownames(flat) <- apply(expand.grid(k=c('c','o'),n=1:numSets), 1,
                        paste0, collapse="")

tdir <- paste0(tempdir(), "/")
write.table(flat, file=paste0(tdir, "testCols.csv"),
            quote=FALSE, row.names = TRUE, col.names=TRUE)

model3 <- mxModel(
  model1,
  mxData(jointdata, 'raw'),
  mxComputeLoop(list(
    LD=mxComputeLoadData(
      'loadData', column=paste0('z',1:2),
      skip.rows=1, skip.cols=1,
      row.names=1, method="oops",
      path=paste0(tdir, "testCols.csv"), verbose=0L),
    mxComputeSetOriginalStarts(),
    mxComputeGradientDescent(),
    mxComputeStandardError(),
    CPT=mxComputeCheckpoint(toReturn=TRUE, standardErrors = TRUE)
  ), i=1:numSets))

expect_error(mxRun(model3), "unknown provider")
model3$compute$steps[['LD']]$method <- 'csv'
model3Fit <- mxRun(model3)

omxCheckEquals(model3Fit$compute$steps[['LD']]$debug$loadCounter, 1L)

discardCols <- c("OpenMxEvals", "iterations", "timestamp",
                 "MxComputeLoop1", "objective", "statusCode", "fitUnits",
                 'testCols.csv:z1', 'testCols.csv:z2')
thr <- c(7, 8, 7, 7, 8, 8, 7, 8, 8, 7, 7, 8, 7, 7, 4, 8, 8, 8, 9, 8,
         9, 9, 11, 10, 10, 9, 9, 8, 10, 10) - 2

log <- model3Fit$compute$steps[['CPT']]$log

omxCheckEquals(log[['testCols.csv:z1']], paste0('c', 1:numSets))
omxCheckEquals(log[['testCols.csv:z2']], paste0('o', 1:numSets))

for (col in discardCols) log[[col]] <- NULL
lmad <- -log10(apply(abs(as.matrix(log - result1)), 2, max))
# names(lmad) <- c()
# cat(deparse(floor(lmad)))
# print(lmad - thr)
omxCheckTrue(all(lmad - thr > 0))

# --------------

totalCol <- ncol(flat) + length(letters)
lcol <- sample.int(totalCol, length(letters))
print(sort(lcol))
flatMap <- (1:totalCol)[-lcol]
flat2 <- matrix("", nrow(flat), totalCol)
flat2[,lcol] <- letters
flat2[,flatMap] <- flat

write.table(flat2, file=paste0(tdir, "testCols.csv"),
            quote=FALSE, row.names = TRUE, col.names=TRUE)

model3$compute$steps$LD$rowFilter <- 1:totalCol %in% lcol
model4Fit <- mxRun(model3)

log <- model4Fit$compute$steps[['CPT']]$log

for (col in discardCols) log[[col]] <- NULL
lmad <- -log10(apply(abs(as.matrix(log - result1)), 2, max))
# names(lmad) <- c()
# cat(deparse(floor(lmad)))
# print(lmad - thr)
omxCheckTrue(all(lmad - thr > 0))
