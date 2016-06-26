# This data is from an email:
#
# Date: Wed, 06 Feb 2013 19:49:24 -0800
# From: Li Cai <lcai at ucla.edu>
# To: Joshua N Pritikin <jpritikin at pobox.com>
# Subject: Re: how did you control item bias in Cai (2010, p. 592) ?

#options(error = browser)
library(OpenMx)
library(rpf)

flexmirt.LL <- 29995.30418

# read data
data.raw <- suppressWarnings(try(read.csv("models/nightly/data/cai2009.csv"), silent=TRUE))
if (is(data.raw, "try-error")) data.raw <- read.csv("data/cai2009.csv")
data.g1 <- as.data.frame(data.raw[data.raw$G==1, 2:13])
data.g2 <- as.data.frame(data.raw[data.raw$G==2, 2:17])

if (0) {
  # for flexmirt
  write.table(data.g1, "cai2009-g1.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(data.g2, "cai2009-g2.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)  
  fm <- read.flexmirt("~/irt/cai2009/cai2009-prm.txt")
  fm$G1$spec <- NULL
  fm$G2$spec <- NULL
}

# from flexMIRT
fm <- structure(list(G1 = structure(list(param = structure(c(0.992675,  0.646717, 0, 0, 0.876469, 1.41764, 1.25402, 0, 0, 0.0826927,  1.76547, 1.20309, 0, 0, -0.346706, 2.1951, 0.844399, 0, 0, -0.978301,  1.37774, 0, 1.06694, 0, 0.992373, 1.80365, 0, 0.814109, 0, 0.213559,  2.15718, 0, 1.58086, 0, -0.418129, 1.18201, 0, 1.56533, 0, -1.24173,  1.80474, 0, 0, 1.0774, 0.810718, 2.60754, 0, 0, 1.23507, 0.0598008,  1.01874, 0, 0, 0.724402, -0.294029, 1.68916, 0, 0, 1.37546, -1.13333 ), .Dim = c(5L, 12L), .Dimnames = list(NULL, c("i1", "i2", "i3",  "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12"))), mean = structure(c(0.822622,  -0.290462, 0.19672, 0.733993), .Names = c("X6", "X7", "X8", "X9" )), cov = structure(c(0.826046, 0, 0, 0, 0, 1.656, 0, 0, 0, 0,  1.11263, 0, 0, 0, 0, 1.07878), .Dim = c(4L, 4L))), .Names = c("param",  "mean", "cov")),
                     G2 = structure(list(param = structure(c(0.992675,  0.646717, 0, 0, 0, 0.876469, 1.41764, 1.25402, 0, 0, 0, 0.0826927,  1.76547, 1.20309, 0, 0, 0, -0.346706, 2.1951, 0.844399, 0, 0,  0, -0.978301, 1.37774, 0, 1.06694, 0, 0, 0.992373, 1.80365, 0,  0.814109, 0, 0, 0.213559, 2.15718, 0, 1.58086, 0, 0, -0.418129,  1.18201, 0, 1.56533, 0, 0, -1.24173, 1.80474, 0, 0, 1.0774, 0,  0.810718, 2.60754, 0, 0, 1.23507, 0, 0.0598008, 1.01874, 0, 0,  0.724402, 0, -0.294029, 1.68916, 0, 0, 1.37546, 0, -1.13333,  1.75531, 0, 0, 0, 1.20652, 0.875564, 1.26308, 0, 0, 0, 1.25013,  0.196607, 1.44526, 0, 0, 0, 0.990354, -0.351181, 1.89461, 0,  0, 0, 0.85611, -1.09382), .Dim = c(6L, 16L), .Dimnames = list(     NULL, c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", "i9",      "i10", "i11", "i12", "i13", "i14", "i15", "i16"))), mean = structure(c(0,  0, 0, 0, 0), .Names = c("X6", "X7", "X8", "X9", "X10")), cov = structure(c(1,  0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 1), .Dim = c(5L, 5L))), .Names = c("param", "mean", "cov" ))), .Names = c("G1", "G2"))

for (col in colnames(data.g1)) data.g1[[col]] <- mxFactor(data.g1[[col]], levels=0:1)
for (col in colnames(data.g2)) data.g2[[col]] <- mxFactor(data.g2[[col]], levels=0:1)

# This function creates a model for a single group.
mk.model <- function(name, data, latent.free) {
  numItems <- dim(data)[2]
  dims <- (1 + numItems/4)
  numPersons <- dim(data)[1]
  spec <- list()
  spec[1:numItems] <- list(rpf.grm(factors = dims))
  
  ip.mat <- mxMatrix(name="item", nrow=dims+1, ncol=numItems,
                     values=c(1.4,rep(0,dims-1),0), free=FALSE)
  rownames(ip.mat) <- c(paste('f', 1:dims, sep=""), "b")
  ip.mat$free[1,] <- TRUE
  ip.mat$free[dims+1,] <- TRUE
  colnames(ip.mat) <- colnames(data)
  for (ix in seq(0,numItems-1,4)) {
    ip.mat$values[2 + ix/4, 1:4 + ix] <- 1
    ip.mat$free[2 + ix/4, 1:4 + ix] <- TRUE
  }
  
  for (ix in 1:numItems) {
    for (px in 1:nrow(ip.mat)) {
      if (!ip.mat$free[px,ix]) next
      ip.mat$labels[px,ix] <- paste(c('p',ix,',',rownames(ip.mat)[px]), collapse='')
    }
  }

  m.mat <- mxMatrix(name="mean", nrow=1, ncol=dims, values=0, free=latent.free)
  colnames(m.mat) <- paste('f', 1:dims, sep="")
  cov.mat.free <- FALSE
  if (latent.free) {
    cov.mat.free <- diag(dims)==1
  }
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims),
                      free=cov.mat.free, lbound=1e-2)
  dimnames(cov.mat) <- list(paste('f', 1:dims, sep=""), paste('f', 1:dims, sep=""))
  
  lname <- paste(name, "latent", sep="")
  latent <- mxModel(lname, m.mat, cov.mat,
		    mxDataDynamic("cov", expectation=paste(name, "expectation", sep=".")),
		    mxExpectationNormal(covariance="cov", means="mean"),
		    mxFitFunctionML())

  m1 <- mxModel(model=name, ip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(
		    ItemSpec=spec,
		    mean=paste(lname, "mean",sep="."),
		    cov=paste(lname, "cov", sep="."),
		    qpoints=21, qwidth=5),
                mxFitFunctionML())
  
  list(ifa=m1, latent=latent)
}

groups <- paste("g", 1:2, sep="")

if (1) {
	# Before fitting the model, check EAP score output against flexMIRT
  g1 <- mk.model("g1", data.g1, TRUE)
  g2 <- mk.model("g2", data.g2, FALSE)
  g1$ifa$item$values[,] <- fm$G1$param
  g1$latent$mean$values <- t(fm$G1$mean)
  g1$latent$cov$values <- fm$G1$cov
  g2$ifa$item$values[,] <- fm$G2$param
  
  # Also check whether we compute the LL correctly given flexMIRT's parameters.
    cModel <- mxModel(model="cModel", c(g1, g2),
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
		      mxComputeOnce('fitfunction', 'fit'))
    cModel.fit <- mxRun(cModel)
    omxCheckCloseEnough(cModel.fit$fitfunction$result, flexmirt.LL, 1e-4)
  
  i1 <- mxModel(cModel,
                mxComputeSequence(steps=list(
                  mxComputeOnce('fitfunction', 'information', "meat"),
                  mxComputeStandardError(),
                  mxComputeHessianQuality())))
  i1 <- mxRun(i1)
  
#  cat(deparse(round(i1$output$standardErrors,3)))
  se <- c(0.085, 0.109, 0.078, 0.131, 0.199, 0.098, 0.148,  0.183, 0.104, 0.165,
          0.134, 0.123, 0.109, 0.149, 0.095, 0.13,  0.123, 0.097, 0.186, 0.23,
          0.124, 0.125, 0.25, 0.138, 0.135,  0.169, 0.101, 0.199, 0.188, 0.127,
          0.084, 0.122, 0.078, 0.146,  0.232, 0.14, 0.104, 0.17, 0.128, 0.174, 0.093,
          0.432, 0.254,  0.324, 0.175, 0.242, 0.125, 0.146, 0.265, 0.1, 0.141, 0.201,
          0.101, 0.189, 0.192, 0.13)
  omxCheckCloseEnough(c(i1$output$standardErrors), se, .01)
  omxCheckCloseEnough(log(i1$output$conditionNumber), 5.1, .2)
}

omxIFAComputePlan <- function(groups) {
  latent.plan <- NULL
  latentFG <- apply(expand.grid(paste(groups,"latent",sep=""), c('mean','cov')), 1, paste, collapse='.')
  if (0) {
    latent.plan <- mxComputeSequence(list(mxComputeOnce(paste(groups, 'expectation', sep='.')),
					  mxComputeOnce(paste(groups, 'expectation', sep='.'),
                                                        "latentDistribution", "copy"),  # c('mean','covariance')
                                          mxComputeOnce('fitfunction', "set-starting-values")),
                                     freeSet=latentFG)
  } else {
	  # default tolerance isn't good enough for stable S-EM results
#    latent.plan <- mxComputeGradientDescent(latentFG, fitfunction="latent.fitfunction")
    latent.plan <- mxComputeNewtonRaphson(latentFG, fitfunction="latent.fitfunction")
  }

  mxComputeSequence(steps=list(
    mxComputeEM(paste(groups, 'expectation', sep='.'), 'scores',
                mxComputeSequence(list(
		    mxComputeNewtonRaphson(freeSet=paste(groups, 'item', sep="."), verbose=0L),
		    latent.plan)),
                #tolerance=1e-10, information="mr1991",
                infoArgs=list(fitfunction=c("fitfunction", "latent.fitfunction")),
		verbose=0L),
    mxComputeStandardError(),
    mxComputeHessianQuality()
  ))
}

latent <- mxModel("latent",
                  mxFitFunctionMultigroup(paste(paste(groups,"latent",sep=""), "fitfunction", sep=".")))

	# Now actually fit the model.
  g1 <- mk.model("g1", data.g1, TRUE)
  g2 <- mk.model("g2", data.g2, FALSE)
  grpModel <- mxModel(model="groupModel", g1, g2, latent,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      omxIFAComputePlan(groups))
  
  #grpModel <- mxOption(grpModel, "Number of Threads", 1)
  
  # NPSOL options:
#   grpModel <- mxOption(grpModel, "Analytic Gradients", 'Yes')
#   grpModel <- mxOption(grpModel, "Verify level", '-1')
#   grpModel <- mxOption(grpModel, "Function precision", '1.0E-7')
  
  grpModel <- mxRun(grpModel)

 omxCheckCloseEnough(grpModel$output$fit, flexmirt.LL, .01)
omxCheckCloseEnough(summary(grpModel)$informationCriteria['AIC:','par'], 30107.30, .01)
omxCheckCloseEnough(summary(grpModel)$informationCriteria['BIC:','par'], 30420.953, .01)

  omxCheckCloseEnough(grpModel$submodels$g2$matrices$item$values,
                      fm$G2$param, .02)
  omxCheckCloseEnough(grpModel$submodels$g1latent$matrices$mean$values, t(fm$G1$mean), .02)
  omxCheckCloseEnough(grpModel$submodels$g1latent$matrices$cov$values, fm$G1$cov, .1)
  
  semse <- c(0.083, 0.104, 0.077, 0.133, 0.198, 0.097, 0.147,  0.191, 0.102, 0.163, 0.131, 0.122,
             0.11, 0.149, 0.093, 0.13,  0.121, 0.095, 0.188, 0.228, 0.122, 0.124, 0.255, 0.139,
             0.137,  0.171, 0.1, 0.207, 0.193, 0.126, 0.086, 0.115, 0.076, 0.142,  0.212, 0.136,
             0.099, 0.156, 0.12, 0.143, 0.096, 0.391, 0.241,  0.284, 0.18, 0.236, 0.126, 0.145,
             0.26, 0.099, 0.142, 0.201,  0.099, 0.178, 0.193, 0.126)
  # cat(deparse(round(grpModel$output$standardErrors,3)))
 print( max(abs(c(grpModel$output$standardErrors) - semse)))
  
  # These are extremely sensitive to small differences in model estimation.
#omxCheckCloseEnough(c(grpModel$output$standardErrors), semse, .1)
#omxCheckCloseEnough(log(grpModel$output$conditionNumber), 5.5, 1)
#omxCheckTrue(grpModel$output$infoDefinite)
  
emstat <- grpModel$compute$steps[[1]]$output
omxCheckCloseEnough(emstat$EMcycles, 57, 16)
#omxCheckCloseEnough(emstat$totalMstep, 334, 40)  # includes latent distribution
#omxCheckCloseEnough(emstat$semProbeCount, 152, 10)

print(grpModel$output$backendTime)

refModels <- mxRefModels(grpModel, run=TRUE)
omxCheckCloseEnough(refModels[['Independence']]$output$fit, 36611.98, .01)
