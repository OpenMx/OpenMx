# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
#library(mvtnorm)

set.seed(7)
correct.LL <- 48939.35

numItems <- 30
numPeople <- 500
maxDim <- 2

items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
	items[[ix]] <- rpf.grm(factors=2)
	correct[[ix]] <- rpf.rparam(items[[ix]])
}
correct.mat <- simplify2array(correct)
correct.mat[1,1:10] <- 0
correct.mat[2,20:30] <- 0

maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))
maxOutcomes <- max(vapply(items, function(i) i@outcomes, 0))

data.g1 <- rpf.sample(numPeople, items, correct.mat)
data.g2 <- rpf.sample(numPeople, items, correct.mat, mean=c(-.5,.8), cov=matrix(c(2,.5,.5,2),nrow=2))
data.g3 <- rpf.sample(numPeople, items, correct.mat, mean=c(-.1,-.8), cov=matrix(c(.9,-.5,-.5,.9),nrow=2))

if (0) {
  # for flexMIRT, write CSV
  write.table(sapply(data.g1, unclass)-1, file="2d-mg-g1.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g2, unclass)-1, file="2d-mg-g2.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g3, unclass)-1, file="2d-mg-g3.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  fm <- read.flexmirt("/home/joshua/irt/ifa-2d-mg/2d-mg-prm.txt")
}

mkgroup <- function(model.name, data, latent.free) {  
  ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                     values=c(1, 1.4, 0), free=TRUE)
  ip.mat@values[correct.mat==0] <- 0
  ip.mat@free[correct.mat==0] <- FALSE
  
  for (ix in 1:numItems) {
    for (px in 1:3) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat@labels[px,ix] <- name
    }
  }
  
  eip.mat <- mxAlgebra(ItemParam, name="EItemParam")
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=0, free=latent.free)
  cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=diag(2), free=latent.free)
  if (latent.free) {
    for (c1 in 1:2) {
      for (c2 in 1:c1) {
        cov.mat@labels[c1,c2] <- paste(model.name, paste(c1, c2, sep=","),sep="-")
        cov.mat@labels[c2,c1] <- paste(model.name, paste(c1, c2, sep=","),sep="-")
      }
    }  
  }
  
  m1 <- mxModel(model=model.name,
                ip.mat, m.mat, cov.mat, eip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=items,
                                  EItemParam="EItemParam", ItemParam="ItemParam",
                                  qpoints=21, qwidth=5),
                mxFitFunctionBA81())
  m1
}

groups <- paste("g", 1:3, sep="")

if (1) {
  fm <- structure(list(G1 = structure(list(param = structure(c(0, 0.547418,  -0.746423, 0, 0.675093, -1.01778, 0, 1.00217, 0.066893, 0, 1.21546,  2.68894, 0, 1.08965, 1.93943, 0, 0.73088, -0.235478, 0, 1.69934,  0.816274, 0, 1.99294, -1.34435, 0, 1.13576, 0.805195, 0, 0.500186,  -0.281224, 0.717699, 1.43936, 0.172262, 0.978424, 0.756688, -0.643423,  1.52223, 0.709748, -0.239081, 1.14345, 2.11652, -0.719694, 0.82372,  0.57575, -0.485341, 0.820995, 2.08929, 0.575516, 0.954974, 1.58002,  -0.28948, 0.837872, 1.13202, 1.51745, 1.27315, 1.32316, -1.56062,  1.20298, 0, -1.05121, 0.971849, 0, 0.49741, 1.31535, 0, 1.01287,  1.63251, 0, 1.24117, 1.48944, 0, -0.943851, 0.664333, 0, -1.71981,  0.46719, 0, -0.518989, 0.848514, 0, 0.371468, 0.924849, 0, -0.365851,  0.641203, 0, -0.303372, 0.884883, 0, 0.600367), .Dim = c(3L,  30L), .Dimnames = list(NULL, c("i1", "i2", "i3", "i4", "i5",  "i6", "i7", "i8", "i9", "i10", "i11", "i12", "i13", "i14", "i15",  "i16", "i17", "i18", "i19", "i20", "i21", "i22", "i23", "i24",  "i25", "i26", "i27", "i28", "i29", "i30"))), mean = structure(c(0,  0), .Names = c("X6", "X7")), cov = structure(c(1, 0, 0, 1), .Dim = c(2L,  2L))), .Names = c("param", "mean", "cov")),
                       G2 = structure(list(param = structure(c(0, 0.547418, -0.746423, 0, 0.675093,      -1.01778, 0, 1.00217, 0.066893, 0, 1.21546, 2.68894, 0, 1.08965,      1.93943, 0, 0.73088, -0.235478, 0, 1.69934, 0.816274, 0,      1.99294, -1.34435, 0, 1.13576, 0.805195, 0, 0.500186, -0.281224,      0.717699, 1.43936, 0.172262, 0.978424, 0.756688, -0.643423,      1.52223, 0.709748, -0.239081, 1.14345, 2.11652, -0.719694,      0.82372, 0.57575, -0.485341, 0.820995, 2.08929, 0.575516,      0.954974, 1.58002, -0.28948, 0.837872, 1.13202, 1.51745,      1.27315, 1.32316, -1.56062, 1.20298, 0, -1.05121, 0.971849,      0, 0.49741, 1.31535, 0, 1.01287, 1.63251, 0, 1.24117, 1.48944,      0, -0.943851, 0.664333, 0, -1.71981, 0.46719, 0, -0.518989,      0.848514, 0, 0.371468, 0.924849, 0, -0.365851, 0.641203,      0, -0.303372, 0.884883, 0, 0.600367), .Dim = c(3L, 30L), .Dimnames = list(         NULL, c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8",          "i9", "i10", "i11", "i12", "i13", "i14", "i15", "i16",          "i17", "i18", "i19", "i20", "i21", "i22", "i23", "i24",          "i25", "i26", "i27", "i28", "i29", "i30"))), mean = structure(c(-0.507038,      0.703947), .Names = c("X6", "X7")), cov = structure(c(1.87718,      0.491825, 0.491825, 2.05385), .Dim = c(2L, 2L))), .Names = c("param",  "mean", "cov")),
                       G3 = structure(list(param = structure(c(0, 0.547418,  -0.746423, 0, 0.675093, -1.01778, 0, 1.00217, 0.066893, 0, 1.21546,  2.68894, 0, 1.08965, 1.93943, 0, 0.73088, -0.235478, 0, 1.69934,  0.816274, 0, 1.99294, -1.34435, 0, 1.13576, 0.805195, 0, 0.500186,  -0.281224, 0.717699, 1.43936, 0.172262, 0.978424, 0.756688, -0.643423,  1.52223, 0.709748, -0.239081, 1.14345, 2.11652, -0.719694, 0.82372,  0.57575, -0.485341, 0.820995, 2.08929, 0.575516, 0.954974, 1.58002,  -0.28948, 0.837872, 1.13202, 1.51745, 1.27315, 1.32316, -1.56062,  1.20298, 0, -1.05121, 0.971849, 0, 0.49741, 1.31535, 0, 1.01287,  1.63251, 0, 1.24117, 1.48944, 0, -0.943851, 0.664333, 0, -1.71981,  0.46719, 0, -0.518989, 0.848514, 0, 0.371468, 0.924849, 0, -0.365851,  0.641203, 0, -0.303372, 0.884883, 0, 0.600367), .Dim = c(3L,  30L), .Dimnames = list(NULL, c("i1", "i2", "i3", "i4", "i5",  "i6", "i7", "i8", "i9", "i10", "i11", "i12", "i13", "i14", "i15",  "i16", "i17", "i18", "i19", "i20", "i21", "i22", "i23", "i24",  "i25", "i26", "i27", "i28", "i29", "i30"))), mean = structure(c(-0.0282586,  -0.827363), .Names = c("X6", "X7")), cov = structure(c(0.891858,  -0.422099, -0.422099, 0.835799), .Dim = c(2L, 2L))), .Names = c("param",  "mean", "cov"))), .Names = c("G1", "G2", "G3"))
  
  g1 <- mkgroup("g1", data.g1, FALSE)
  g2 <- mkgroup("g2", data.g2, TRUE)
  g3 <- mkgroup("g3", data.g3, TRUE)
  
  g1@matrices$ItemParam@values <- fm$G1$param
  g2@matrices$ItemParam@values <- fm$G2$param
  g3@matrices$ItemParam@values <- fm$G3$param
  g2@matrices$mean@values <- t(fm$G2$mean)
  g3@matrices$mean@values <- t(fm$G3$mean)
  g2@matrices$cov@values <- fm$G2$cov
  g3@matrices$cov@values <- fm$G3$cov
  
  cModel <- mxModel(model="cModel", g1,g2,g3,
                    mxComputeOnce(paste(groups, 'expectation', sep='.'), context='EM'))
  for (grp in groups) cModel@submodels[[grp]]@expectation@scores <- 'full'
  cModel.eap <- mxRun(cModel)
  
  fm.sco <- suppressWarnings(try(read.table("models/nightly/data/2d-mg-sco.txt"), silent=TRUE))
  if (is(fm.sco, "try-error")) fm.sco <- read.table("data/2d-mg-sco.txt")
  colnames(fm.sco) <- c("group", "row", "s1", "s2", "se1", "se2", paste("cov",1:3,sep=""))
  
  omxCheckCloseEnough(as.matrix(fm.sco[fm.sco$group==1,-1:-2]),
                      cModel.eap@submodels$g1@expectation@scores.out, 1e-3)
  omxCheckCloseEnough(as.matrix(fm.sco[fm.sco$group==2,-1:-2]),
                      cModel.eap@submodels$g2@expectation@scores.out, 1e-3)
  omxCheckCloseEnough(as.matrix(fm.sco[fm.sco$group==3,-1:-2]),
                      cModel.eap@submodels$g3@expectation@scores.out, 1e-3)

  cModel <- mxModel(cModel,
                    mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                    mxComputeSequence(steps=list(
                      mxComputeOnce(paste(groups, 'expectation', sep='.')),
                      mxComputeOnce('fitfunction'))))
  for (grp in groups) cModel@submodels[[grp]]@expectation@scores <- 'omit'
  cModel.fit <- mxRun(cModel)
  omxCheckCloseEnough(cModel.fit@output$minimum, correct.LL, .01)
}

if (1) {
  g1 <- mkgroup("g1", data.g1, FALSE)
  g2 <- mkgroup("g2", data.g2, TRUE)
  g3 <- mkgroup("g3", data.g3, TRUE)
  
  grpModel <- mxModel(model="groupModel", g1, g2, g3,
                      mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                      mxComputeIterate(steps=list(
                        mxComputeOnce(paste(groups, "EItemParam", sep=".")),
                        mxComputeOnce(paste(groups, 'expectation', sep='.'), context='EM'),
                        mxComputeNewtonRaphson(free.set=paste(groups, 'ItemParam', sep=".")),
                        mxComputeOnce(paste(groups, 'expectation', sep=".")),
#                        mxComputeOnce('fitfunction', start=TRUE,
#                                      free.set=apply(expand.grid(groups, c('mean','cov')), 1, paste, collapse='.'))
			mxComputeGradientDescent(start=TRUE, useGradient=TRUE,
						 free.set=apply(expand.grid(groups, c('mean','cov')), 1, paste, collapse='.'))
                      ), verbose=TRUE))

  grpModel <- mxOption(grpModel, "Analytic Gradients", 'Yes')
	grpModel <- mxOption(grpModel, "Verify level", '-1')
  grpModel <- mxOption(grpModel, "Function precision", '1.0E-5')

  grpModel <- mxRun(grpModel, silent=TRUE)
  omxCheckCloseEnough(grpModel@output$minimum, correct.LL, .01)
  omxCheckCloseEnough(c(grpModel@submodels$g2@matrices$mean@values), c(-.507, .703), .01)
  omxCheckCloseEnough(c(grpModel@submodels$g2@matrices$cov@values), c(1.877, .491, .491, 2.05), .01)
  omxCheckCloseEnough(c(grpModel@submodels$g3@matrices$mean@values), c(-.028, -.827), .01)
  omxCheckCloseEnough(c(grpModel@submodels$g3@matrices$cov@values), c(.892, -.422, -.422, .836), .01)
}
