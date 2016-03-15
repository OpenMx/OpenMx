library(OpenMx)

numManifestsPerLatent <- 5

mkModel <- function(shuffle, fellner) {
  varStruct <- expand.grid(l=1:3, i=1:numManifestsPerLatent)
  manifestVars <- paste0('i', apply(varStruct, 1, paste0, collapse=''))
  if (shuffle) {
    manifestVars <- sample(manifestVars, length(manifestVars))
  }
  latentVars <- paste0('l', unique(varStruct$l))
  allVars <- c(manifestVars, latentVars)
  if (shuffle) {
    allVars <- sample(allVars, length(allVars))
  }
  
  Fval <- diag(length(manifestVars))[,match(allVars, manifestVars)]
  Fval[is.na(Fval)] <- 0
  
  ta1 <- mxModel(
    model="tangle",
    mxMatrix("Full", length(manifestVars), length(allVars),
             values=Fval,
             dimnames=list(manifestVars, allVars), name="F"),
    mxMatrix("Symm", length(allVars), length(allVars),
             values=diag(length(allVars)),
             free=diag(length(allVars)) == 1,
             dimnames=list(allVars, allVars), name="S"),
    mxMatrix("Full", length(allVars), length(allVars),
             values=0,
             dimnames=list(allVars, allVars), name="A"),
    mxMatrix("Full", 1, length(allVars),
             free=!is.na(match(allVars, manifestVars)),
             dimnames=list(NULL, allVars), name="M"),
    mxExpectationRAM(M="M"),
    mxFitFunctionML(fellner=fellner))
  
  for (lx in 1:length(latentVars)) {
    lvar <- paste0('l', lx)
    ivar <- paste0(paste0('i', lx), 1:numManifestsPerLatent)
    ta1$A$values[ivar, lvar] <- 1
    ta1$A$free[ivar, lvar[-1]] <- TRUE
  }
  
  ta1$S$free[latentVars, latentVars] <- TRUE
  ta1
}

set.seed(1)
ta1 <- mxGenerateData(mkModel(FALSE, FALSE), nrow=25, returnModel=TRUE)

fit1 <- mxRun(mxModel(mkModel(TRUE, TRUE), ta1$data))
fit2 <- mxRun(mxModel(mkModel(FALSE, TRUE), ta1$data))

omxCheckCloseEnough(fit1$output$fit - fit2$output$fit, 0, 1e-8)

omxCheckEquals(names(fit1$expectation$debug$g1$fullMean[1:ncol(fit1$F)]),
               paste0('tangle.', colnames(fit1$F)))

omxCheckTrue(length(rle(fit1$expectation$debug$g1$latentFilter)$lengths) !=
               length(rle(fit2$expectation$debug$g1$latentFilter)$lengths))

if (0) {
  # F does two things:
  # drops out latent variables
  # permutes the manifest variables
  
  all(rownames(fit1$F$values) == fit1$expectation$.runDims)
  filter <- !apply(fit1$F$values, 2, function(x) all(x==0))
  
  map1 <- apply(fit1$F$values != 0, 2, function(x) ifelse(any(x), which(x), NA))
  man1 <- colnames(fit1$F)[filter]
  # permit F input order to F output order
  all(man1[ match(rownames(fit1$F), man1) ] == rownames(fit1$F))
  # permute data to F input order
  all(rownames(fit1$F)[ match(man1, rownames(fit1$F)) ] == man1)
  all(match(man1, rownames(fit1$F)) == map1[!is.na(map1)])
}
