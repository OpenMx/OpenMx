# This data is from an email:
#
# Date: Wed, 06 Feb 2013 19:49:24 -0800
# From: Li Cai <lcai@ucla.edu>
# To: Joshua N Pritikin <jpritikin@pobox.com>
# Subject: Re: how did you control item bias in Cai (2010, p. 592) ?

#options(error = utils::recover)
library(OpenMx)
library(rpf)

correct.LL <- 29995.30418  # from flexMIRT

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
mk.model <- function(model.name, data, latent.free) {
  numItems <- dim(data)[2]
  numPersons <- dim(data)[1]
  spec <- list()
  spec[1:numItems] <- rpf.grm(factors = 2)
  
  dims <- (1 + numItems/4)
  design <- matrix(c(rep(1L,numItems),
                     as.integer(kronecker(2:dims,rep(1,4)))), byrow=TRUE, ncol=numItems)
  
  ip.mat <- mxMatrix(name="ItemParam", nrow=3, ncol=numItems,
                     values=c(1.4,1,0),
                     free=c(TRUE,TRUE,TRUE))
  
  for (ix in 1:numItems) {
    for (px in 1:3) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat@labels[px,ix] <- name
    }
  }

  m.mat <- mxMatrix(name="mean", nrow=1, ncol=dims, values=0, free=latent.free)
  cov.mat.free <- FALSE
  if (latent.free) {
    cov.mat.free <- diag(dims)==1
  }
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims),
                      free=cov.mat.free)
  
  m1 <- mxModel(model=model.name, ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(
                  ItemSpec=spec,
                  design=design,
                  ItemParam="ItemParam",
                  mean="mean", cov="cov",
                  qpoints=21, qwidth=5),
                mxFitFunctionML())
  m1
}

groups <- paste("g", 1:2, sep="")

if (1) {
	# Before fitting the model, check EAP score output against flexMIRT
  g1 <- mk.model("g1", data.g1, TRUE)
  g2 <- mk.model("g2", data.g2, FALSE)
  g1@matrices$ItemParam@values <-
    rbind(fm$G1$param[1,], apply(fm$G1$param[2:4,], 2, sum), fm$G1$param[5,])
  g1@matrices$mean@values <- t(fm$G1$mean)
  g1@matrices$cov@values <- fm$G1$cov
  g2@matrices$ItemParam@values <-
    rbind(fm$G2$param[1,], apply(fm$G2$param[2:5,], 2, sum), fm$G2$param[6,])
  
  cModel <- mxModel(model="cModel", g1,g2,
                    mxComputeOnce(paste(groups, 'expectation', sep='.'), context='EM'))
#  cModel <- mxOption(cModel, "Number of Threads", 1)
  for (grp in groups) cModel@submodels[[grp]]@expectation@scores <- 'full'
  cModel.eap <- mxRun(cModel)

  fm.sco.g1 <- suppressWarnings(try(read.table("models/nightly/data/cai2009-sco-g1.txt"), silent=TRUE))
  if (is(fm.sco.g1, "try-error")) fm.sco.g1 <- read.table("data/cai2009-sco-g1.txt")
  fm.sco.g2 <- suppressWarnings(try(read.table("models/nightly/data/cai2009-sco-g2.txt"), silent=TRUE))
  if (is(fm.sco.g2, "try-error")) fm.sco.g2 <- read.table("data/cai2009-sco-g2.txt")
  colnames(fm.sco.g1) <- c("grp","id",colnames(cModel.eap@submodels$g1@expectation@output$scores))
  colnames(fm.sco.g2) <- c("grp","id",colnames(cModel.eap@submodels$g2@expectation@output$scores))
  
  scores.g1 <- cModel.eap@submodels$g1@expectation@output$scores
  omxCheckCloseEnough(as.matrix(fm.sco.g1[,-1:-2]),
                      scores.g1, 1e-3)
  omxCheckCloseEnough(as.matrix(fm.sco.g2[,-1:-2]),
                      cModel.eap@submodels$g2@expectation@output$scores, 1e-3)

  # Also check whether we compute the LL correctly given flexMIRT's parameters.
    cModel <- mxModel(cModel,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      mxComputeSequence(steps=list(
                        mxComputeOnce(paste(groups, 'expectation', sep=".")),
                        mxComputeOnce('fitfunction', fit=TRUE,
				      free.set=apply(expand.grid(groups, c('mean','cov')), 1, paste, collapse='.')))))
    cModel.fit <- mxRun(cModel)
    omxCheckCloseEnough(cModel.fit@fitfunction@result, correct.LL, 1e-4)
  
  i1 <- mxModel(cModel,
                mxComputeSequence(steps=list(
                  mxComputeOnce(paste(groups, 'expectation', sep='.')),
                  mxComputeOnce('fitfunction', information=TRUE, info.method="meat"),
                  mxComputeStandardError(),
                  mxComputeHessianQuality())))
  i1 <- mxRun(i1)
  
#  cat(deparse(round(i1@output$standardErrors,3)))
  se <- c(0.085, 0.109, 0.078, 0.131, 0.199, 0.098, 0.148,  0.183, 0.104, 0.165,
          0.134, 0.123, 0.109, 0.149, 0.095, 0.13,  0.123, 0.097, 0.186, 0.23,
          0.124, 0.125, 0.25, 0.138, 0.135,  0.169, 0.101, 0.199, 0.188, 0.127,
          0.084, 0.122, 0.078, 0.146,  0.232, 0.14, 0.104, 0.17, 0.128, 0.174, 0.093,
          0.432, 0.254,  0.324, 0.175, 0.242, 0.125, 0.146, 0.265, 0.1, 0.141, 0.201,
          0.101, 0.189, 0.192, 0.13)
  omxCheckCloseEnough(c(i1@output$standardErrors), se, .01)
  omxCheckCloseEnough(i1@output$conditionNumber, 199, 1) 
}

omxIFAComputePlan <- function(groups) {
  mxComputeSequence(steps=list(
    mxComputeEM(paste(groups, 'expectation', sep='.'),
                mxComputeNewtonRaphson(free.set=paste(groups, 'ItemParam', sep=".")),
                mxComputeOnce('fitfunction', fit=TRUE,
                              free.set=apply(expand.grid(groups, c('mean','cov')), 1, paste, collapse='.')),
                tolerance=1e-5, information=TRUE),
    mxComputeStandardError(),
    mxComputeHessianQuality()
  ))
}

	# Now actually fit the model.
  g1 <- mk.model("g1", data.g1, TRUE)
  g2 <- mk.model("g2", data.g2, FALSE)
  grpModel <- mxModel(model="groupModel", g1, g2,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      omxIFAComputePlan(groups))
  
  #grpModel <- mxOption(grpModel, "Number of Threads", 1)
  
  # NPSOL options:
#   grpModel <- mxOption(grpModel, "Analytic Gradients", 'Yes')
#   grpModel <- mxOption(grpModel, "Verify level", '-1')
#   grpModel <- mxOption(grpModel, "Function precision", '1.0E-7')
  
  grpModel <- mxRun(grpModel)
    
  omxCheckCloseEnough(grpModel@output$minimum, correct.LL, .01)
  omxCheckCloseEnough(grpModel@submodels$g2@matrices$ItemParam@values,
                      rbind(fm$G2$param[1,], apply(fm$G2$param[2:5,], 2, sum), fm$G2$param[6,]), .01)
  omxCheckCloseEnough(grpModel@submodels$g1@matrices$mean@values, t(fm$G1$mean), .01)
  omxCheckCloseEnough(grpModel@submodels$g1@matrices$cov@values, fm$G1$cov, .01)
  
  semse <- c(0.084, 0.105, 0.078, 0.134, 0.201, 0.1, 0.154, 0.196,  0.106, 0.176,
             0.132, 0.121, 0.11, 0.149, 0.094, 0.133, 0.122,  0.092, 0.199, 0.232,
             0.125, 0.125, 0.26, 0.141, 0.138, 0.171,  0.102, 0.214, 0.193, 0.128,
             0.086, 0.115, 0.078, 0.142, 0.211,  0.14, 0.102, 0.162, 0.127, 0.156,
             0.1, 0.408, 0.243, 0.289, 0.181,  0.238, 0.128, 0.146, 0.262, 0.1, 0.143,
             0.202, 0.1, 0.181, 0.193,  0.129)
  # cat(deparse(round(grpModel@output$standardErrors,3)))
  # max(abs(c(grpModel@output$standardErrors) - semse))
  
  # These are extremely sensitive to small differences in model estimation.
  omxCheckCloseEnough(c(grpModel@output$standardErrors), semse, .05)
  omxCheckCloseEnough(grpModel@output$conditionNumber, 250, 100)
  
  print(grpModel@output$backendTime)
