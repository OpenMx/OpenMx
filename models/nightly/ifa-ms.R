library(OpenMx)
library(rpf)

set.seed(1)
m2.data <- suppressWarnings(try(read.table("models/nightly/data/ms-data.csv"), silent=TRUE))
if (is(m2.data, "try-error")) m2.data <- read.table("data/ms-data.csv")
m2.data[m2.data==-9] <- NA
m2.data <- m2.data + 1

m2.data <- data.frame(lapply(m2.data, mxFactor, levels=1:5))

gpcm <- function(outcomes) {
  rpf.nrm(outcomes, T.c=lower.tri(diag(outcomes-1),TRUE) * -1)
  #   rpf.nrm(outcomes, T.c=diag(outcomes-1))
}

m2.spec <- list()
m2.spec[1:22] <- gpcm(5)
m2.spec[2] <- gpcm(4)
m2.spec[5] <- gpcm(3)
m2.spec[6] <- gpcm(4)
m2.spec[13:14] <- gpcm(4)

m2.numItems <- length(m2.spec)
m2.maxParam <-max(sapply(m2.spec, rpf.numParam))

ip.mat <- mxMatrix(name="ItemParam", nrow=m2.maxParam, ncol=m2.numItems,
                   values=c(1, 1, rep(0, m2.maxParam-2)), free=FALSE)
ip.mat@labels[1,] <- 'a1'
ip.mat@free[1,] <- TRUE
rstart <- lapply(m2.spec, rpf.rparam)
for (ix in 1:m2.numItems) {
  thr <- m2.spec[[ix]]@outcomes - 1
  ip.mat@free[(2+thr):(1+2*thr), ix] <- TRUE
  ip.mat@values[ip.mat@free[,ix],ix] <- rstart[[ix]][ip.mat@free[,ix]]
}
ip.mat@values[!is.na(ip.mat@labels) & ip.mat@labels == 'a1'] <-
  sample(ip.mat@values[!is.na(ip.mat@labels) & ip.mat@labels == 'a1'], 1)

#  m2.fmfit <- read.flexmirt("~/2012/sy/fm/ms-rasch-prm.txt")
# cat(deparse(round(m2.fmfit$G1$param,6)))
fmfit <- structure(c(0.941583, 1, 0, 0, 0, -0.676556, 0.758794, -0.802595,  1.28891, 0.941583, 1, 0, 0, -0.182632, 0.897435, 1.30626, NA,  NA, 0.941583, 1, 0, 0, 0, 0.177835, -1.82185, 0.005832, -0.81109,  0.941583, 1, 0, 0, 0, -1.15962, -1.229, 0.032677, 0.4922, 0.941583,  1, 0, 0.457533, 0.324595, NA, NA, NA, NA, 0.941583, 1, 0, 0,  -2.69186, -1.04012, 1.61232, NA, NA, 0.941583, 1, 0, 0, 0, -1.38231,  0.034368, -1.214, -0.648291, 0.941583, 1, 0, 0, 0, -1.85655,  -1.17135, -0.262079, -0.531158, 0.941583, 1, 0, 0, 0, -1.29475,  -0.376539, 0.02024, 0.135187, 0.941583, 1, 0, 0, 0, -1.38279,  0.524151, -0.508742, 0.633671, 0.941583, 1, 0, 0, 0, -0.979595,  -0.048528, 0.659669, 0.544857, 0.941583, 1, 0, 0, 0, -2.09039,  -1.45472, -0.472137, -0.666386, 0.941583, 1, 0, 0, 0.174682,  0.645437, 0.907132, NA, NA, 0.941583, 1, 0, 0, -0.842216, 0.490717,  1.28034, NA, NA, 0.941583, 1, 0, 0, 0, -0.913355, -0.319602,  -0.310164, -0.15536, 0.941583, 1, 0, 0, 0, 0.567085, -1.56762,  0.884553, 0.122113, 0.941583, 1, 0, 0, 0, -0.152985, -0.341317,  -0.183837, 1.17952, 0.941583, 1, 0, 0, 0, 0.168869, -0.490354,  0.373892, 1.29714, 0.941583, 1, 0, 0, 0, -0.827385, 0.626197,  -1.52994, 0.494209, 0.941583, 1, 0, 0, 0, 0.511263, -0.750358,  1.01852, 0.840026, 0.941583, 1, 0, 0, 0, 0.968905, -0.009671,  1.52297, 1.69255, 0.941583, 1, 0, 0, 0, 1.89582, 0.051828, 2.25758,  1.52469), .Dim = c(9L, 22L), .Dimnames = list(NULL, c("i1", "i2",  "i3", "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12",  "i13", "i14", "i15", "i16", "i17", "i18", "i19", "i20", "i21",  "i22")))
#  ip.mat@values <- m2.fmfit$G1$param

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

if (1) {
  cip.mat <- ip.mat
  cip.mat@values <- fmfit
  cM <- mxModel(model="ms", m.mat, cov.mat, cip.mat,
                mxData(observed=m2.data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=m2.spec,
                                  ItemParam="ItemParam"),
                mxFitFunctionML(),
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation'),
                  mxComputeOnce('fitfunction', fit=TRUE, free.set=c("mean", "cov"))
                )))
  cM <- mxRun(cM, silent=TRUE)
  omxCheckCloseEnough(cM@fitfunction@result, 50661.38, .01)
}

m2 <- mxModel(model="m2", m.mat, cov.mat, ip.mat,
              mxData(observed=m2.data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec=m2.spec,
                                ItemParam="ItemParam"),
              mxFitFunctionML(),
              mxComputeIterate(steps=list(
                mxComputeOnce('expectation', context='EM'),
#                mxComputeGradientDescent(free.set='ItemParam', useGradient=TRUE),
                mxComputeNewtonRaphson(free.set='ItemParam'),
                mxComputeOnce('expectation'),
                mxComputeOnce('fitfunction', fit=TRUE, free.set=c("mean", "cov"))
              )))
#  m2 <- mxOption(m2, "Number of Threads", 1)
m2 <- mxRun(m2, silent=TRUE)
omxCheckCloseEnough(m2@fitfunction@result, 50661.377, .01)

#print(m2@matrices$ItemParam@values - fmfit)
print(m2@output$backendTime)
