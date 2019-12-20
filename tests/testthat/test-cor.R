library(OpenMx)
library(testthat)

data(demoOneFactor)

myData   <- mxData(cor(demoOneFactor), type = "cor", numObs = 500)
manifests = names(demoOneFactor)
latents <- c("G")

paths <- list(mxPath(from=latents, to=manifests),
              mxPath(from=manifests, arrows=2),
              mxPath(from=latents, arrows=2, free=FALSE, values=1.0))

fl <- mxModel("L", 
              type="LISREL",
              manifestVars=list(exo=manifests), 
              latentVars=list(exo=latents),
              paths, myData)

fr <- mxModel("R",
              type="RAM",
              manifestVars = manifests,
              latentVars = latents,
              paths, myData)

mg <- mxModel("Both", fl, fr,
        mxFitFunctionMultigroup(c('L','R')))

fl1 <- mxRun(fl)
fr1 <- mxRun(fr)
mg1 <- mxRun(mg)

for (fit in list(fl1, fr1, mg1$L, mg1$R)) {
  expect_equal(diag(fit$expectedCovariance$values), rep(1,5))
  expect_equivalent(diag(mxGetExpected(fit, 'covariance')), rep(1,5))
}

expect_equivalent(fl1$output$constraintJacobian[,paste0("L.TD[",1:5,",",1:5,"]")], diag(5))
expect_equivalent(fr1$output$constraintJacobian[,paste0("R.S[",1:5,",",1:5,"]")], diag(5))

expect_equivalent(mg1$output$constraintJacobian,
                  rbind(
                    cbind(fl1$output$constraintJacobian, matrix(0, 5,10)),
                    cbind(matrix(0,5,10), fr1$output$constraintJacobian)),
                  tolerance=1e-6)

fr1$expectedCovariance$free[1,1] <- TRUE
expect_error(mxRun(fr1), "Free parameters are not allowed")
fr1$expectedCovariance$free[1,1] <- FALSE

fr1$expectedCovariance$labels[1,1] <- "bob"
expect_error(mxRun(fr1), "Labels are not allowed")
fr1$expectedCovariance$labels[1,1] <- NA

expect_error(mxRun(mxModel(fr,
                           mxMatrix(name="expectedCovariance",nrow=1,ncol=1))),
             "an object named 'expectedCovariance' already exists")

# ---- #

fm <- mxModel("One Factor", type="RAM",
              manifestVars = manifests,
              latentVars = latents, paths,
              mxPath(from = 'one', to = manifests),
              mxData(demoOneFactor, type = "raw"), 
              mxMatrix(nrow=1, ncol=1, name="expectedMean"))
fm$expectation$expectedMean <- "expectedMean"
expect_error(mxRun(fm), "Matrix 'expectedMean' must be dimension 1x5")

fm <- mxModel(fm, mxMatrix(nrow=1, ncol=5, name="expectedMean"))
fm <- mxRun(fm)
expect_equivalent(fm$expectedMean$values, colMeans(demoOneFactor), tolerance=1e-5)

# do constructor args work?
fm <- mxModel(
  fm,
  mxMatrix(nrow=5, ncol=5, name="expectedCovariance"),
  mxExpectationRAM(M="M", expectedMean="expectedMean",
                   expectedCovariance="expectedCovariance"))
fm <- mxRun(fm)
expect_equivalent(fm$expectedMean$values, colMeans(demoOneFactor), tolerance=1e-5)
expect_equivalent(max(abs(fm$expectedCovariance$values - cov(demoOneFactor))),
                  0, tolerance=1e-2)
