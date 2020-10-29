library(OpenMx)

set.seed(1)

data("jointdata", package ="OpenMx", verbose= TRUE)
jointData <- jointdata

jointData[,c(2,4,5)] <- mxFactor(jointData[,c(2,4,5)], 
                                 levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

jointData$z1c <- with(jointData, z1 * .1 + rnorm(nrow(jointData)))

jointData$z2c <- with(jointData, rnorm(nrow(jointData), mean=unclass(z2)*.2))

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
colnames(thresh) <- paste0('z', c(2,4,5))

unweighted <- rbind(
  jointData[1,],
  jointData[1,],
  jointData[1,],
  jointData[1,],
  jointData[-2:-3,])

weighted <- cbind(jointData, freq=c(5L,0L,0L,rep(1L,nrow(jointData)-3)))

#------------------------------------------------------------------------------
# Model definition

uwModel <- mxModel(
  "JointRAM", type="RAM", thresh,
  manifestVars = paste0('z', 1:5),
  latentVars = c('G','z1c','z2c'),
  mxData(unweighted, 'raw'),
  mxPath('one', paste0('z', c(1,3)), free=TRUE),
  mxPath(paste0('z', c(1,3)), arrows=2, free=TRUE, values=.5),
  mxPath(paste0('z', c(2,4,5)), arrows=2, free=FALSE, values=.5),
  mxPath('G', arrows=2, values=1, free=FALSE),
  mxPath('G', paste0('z', 1:5), free=TRUE, values=1, lbound=0),
  mxPath('one', 'z1c', free=FALSE, labels="data.z1c"),
  mxPath('one', 'z2c', free=FALSE, labels="data.z2c"),
  mxPath('z1c', 'z1'),
  mxPath('z2c', 'z2'),
  mxFitFunctionWLS())

uwModel$expectation$thresholds <- 'T'

wModel <- mxModel(uwModel, mxData(weighted, 'raw', verbose=0L, frequency="freq"))

uwModel <- mxRun(uwModel)
wModel <- mxRun(wModel)

os1 <- uwModel$data$observedStats
os2 <- wModel$data$observedStats

omxCheckCloseEnough(vech(os1$cov), vech(os2$cov), 1e-12)
omxCheckCloseEnough(os1$means, os2$means, 1e-12)
omxCheckCloseEnough(os1$slope, os2$slope, 1e-12)
mask <- !is.na(os1$thresholds)
omxCheckCloseEnough(os1$thresholds[mask], os2$thresholds[mask], 1e-12)
omxCheckCloseEnough(vech(os1$acov), vech(os2$acov), 1e-10)

#omxCheckCloseEnough(cor(coef(uwModel), coef(wModel)), 1, 1e-9)

omxCheckCloseEnough(uwModel$output$standardErrors,
                    wModel$output$standardErrors, 1e-3)

wald.se <- c(uwModel$output$standardErrors)

# -----------

contTestData <- Bollen[,1:8]
unweighted <- rbind(
  contTestData[1,],
  contTestData[1,],
  contTestData[1,],
  contTestData[1,],
  contTestData[-2:-3,])

weighted <- cbind(contTestData, freq=c(5L,0L,0L,rep(1L,nrow(contTestData)-3)))

uwData <- omxAugmentDataWithWLSSummary(mxData(unweighted, 'raw'))
wData <- omxAugmentDataWithWLSSummary(mxData(weighted, 'raw', frequency="freq"))

os1 <- uwData$observedStats
os2 <- wData$observedStats
omxCheckCloseEnough(os1$cov, os2$cov, 1e-12)
omxCheckCloseEnough(os1$fullWeight, os2$fullWeight, 2e-12)

# -----------

uwModelB <- mxBootstrap(uwModel)
uw.sum <- summary(uwModelB)
omxCheckCloseEnough(cor(uw.sum$bootstrapSE, wald.se), 1, .15)

omxCheckCloseEnough((uw.sum$bootstrapSE -  wald.se) / uw.sum$bootstrapSE, rep(0,length(wald.se)), .41)
