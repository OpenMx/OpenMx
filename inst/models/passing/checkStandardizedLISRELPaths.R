library(OpenMx)

# Set up data
vNames <- paste("v", as.character(1:6), sep="")
dimList <- list(vNames, vNames)
covData <- matrix(
  c(0.9223099, 0.1862938, 0.4374359, 0.8959973, 0.9928430, 0.5320662,
    0.1862938, 0.2889364, 0.3927790, 0.3321639, 0.3371594, 0.4476898,
    0.4374359, 0.3927790, 1.0069552, 0.6918755, 0.7482155, 0.9013952,
    0.8959973, 0.3321639, 0.6918755, 1.8059956, 1.6142005, 0.8040448,
    0.9928430, 0.3371594, 0.7482155, 1.6142005, 1.9223567, 0.8777786,
    0.5320662, 0.4476898, 0.9013952, 0.8040448, 0.8777786, 1.3997558
    ), nrow=6, ncol=6, byrow=TRUE, dimnames=dimList)

# Create LISREL matrices
mLX <- mxMatrix("Full", values=c(.5, .6, .8, rep(0, 6), .4, .7, .5),
          name="LX", nrow=6, ncol=2,
          free=c(TRUE,TRUE,TRUE,rep(FALSE, 6),TRUE,TRUE,TRUE),
          dimnames=list(vNames, c("x1","x2")))
mTD <- mxMatrix("Diag", values=c(rep(.2, 6)), 
          name="TD", nrow=6, ncol=6, free=TRUE,
          dimnames=dimList)
mPH <- mxMatrix("Symm", values=c(1, .3, 1), 
          name="PH", nrow=2, ncol=2, free=c(FALSE, TRUE, FALSE),
          dimnames=list(c("x1","x2"),c("x1","x2")))

expFunction <- mxExpectationLISREL(LX="LX", TD="TD", PH="PH")
tmpData <- mxData(observed=covData, type="cov", numObs=100)
fitFunction <- mxFitFunctionML()
tmpModel <- mxModel(model="exampleModel",
                    mLX, mTD, mPH, expFunction, fitFunction, tmpData)
tmpModelOut <- mxRun(tmpModel)

# Standardize LISREL paths
out <- mxStandardizeLISRELpaths(tmpModelOut, SE=TRUE)

# Manual calculations to verify
LX_val <- tmpModelOut$LX$values
TD_val <- tmpModelOut$TD$values
PH_val <- tmpModelOut$PH$values

cov_X <- LX_val %*% PH_val %*% t(LX_val) + TD_val
sd_X <- sqrt(diag(cov_X))
sd_xi <- sqrt(diag(PH_val))

# Standardize LX
LX_std_manual <- diag(1/sd_X) %*% LX_val %*% diag(sd_xi)
# Standardize TD
TD_std_manual <- diag(1/sd_X) %*% TD_val %*% diag(1/sd_X)
# Standardize PH
PH_std_manual <- diag(1/sd_xi) %*% PH_val %*% diag(1/sd_xi)

# Verify LX elements in output
for (i in 1:nrow(LX_std_manual)) {
  for (j in 1:ncol(LX_std_manual)) {
    val <- LX_std_manual[i, j]
    if (val != 0) {
      row_idx <- which(out$matrix == "LX" & out$row == vNames[i] & out$col == c("x1", "x2")[j])
      omxCheckCloseEnough(out$Std.Value[row_idx], val, 1e-6)
    }
  }
}

# Verify TD elements
for (i in 1:nrow(TD_std_manual)) {
  val <- TD_std_manual[i, i]
  if (val != 0) {
    row_idx <- which(out$matrix == "TD" & out$row == vNames[i] & out$col == vNames[i])
    omxCheckCloseEnough(out$Std.Value[row_idx], val, 1e-6)
  }
}

# Verify PH elements
PH_names <- c("x1", "x2")
for (i in 1:2) {
  for (j in 1:i) {
    val <- PH_std_manual[i, j]
    if (val != 0) {
      row_idx <- which(out$matrix == "PH" & out$row == PH_names[i] & out$col == PH_names[j])
      omxCheckCloseEnough(out$Std.Value[row_idx], val, 1e-6)
    }
  }
}

# Check that standard errors are non-null and numeric
omxCheckTrue(all(!is.na(out$Std.SE)))
omxCheckTrue(all(is.numeric(out$Std.SE)))

print("All LISREL path standardization tests passed!")
