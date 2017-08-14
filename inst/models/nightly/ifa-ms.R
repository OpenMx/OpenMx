library(OpenMx)
library(rpf)

set.seed(1)
m2.data <- suppressWarnings(try(read.table("models/nightly/data/ms-data.csv"), silent=TRUE))
if (is(m2.data, "try-error")) m2.data <- read.table("data/ms-data.csv")
m2.data[m2.data==-9] <- NA
m2.data <- m2.data + 1

gpcm <- function(outcomes) {
  rpf.nrm(outcomes, T.c=lower.tri(diag(outcomes-1),TRUE) * -1)
  #   rpf.nrm(outcomes, T.c=diag(outcomes-1))
}

m2.spec <- list()
m2.spec[1:22] <- list(gpcm(5))
m2.spec[2] <- list(gpcm(4))
m2.spec[5] <- list(gpcm(3))
m2.spec[6] <- list(gpcm(4))
m2.spec[13:14] <- list(gpcm(4))

m2.numItems <- length(m2.spec)

for (c in 1:m2.numItems) {
  m2.data[[c]] <- mxFactor(m2.data[[c]], levels=1:m2.spec[[c]]$outcomes)
}

m2.maxParam <-max(sapply(m2.spec, rpf.numParam))

ip.mat <- mxMatrix(name="item", nrow=m2.maxParam, ncol=m2.numItems,
                   values=c(1, 1, rep(0, m2.maxParam-2)), free=FALSE)
colnames(ip.mat) <- colnames(m2.data)
rownames(ip.mat) <- c('f1', rep('n', nrow(ip.mat)-1))
ip.mat$labels[1,] <- 'a1'
ip.mat$free[1,] <- TRUE
rstart <- lapply(m2.spec, rpf.rparam, version=1)
for (ix in 1:m2.numItems) {
  thr <- m2.spec[[ix]]$outcomes - 1
  ip.mat$free[(2+thr):(1+2*thr), ix] <- TRUE
  ip.mat$values[ip.mat$free[,ix],ix] <- rstart[[ix]][ip.mat$free[,ix]]
}
ip.mat$values[!is.na(ip.mat$labels) & ip.mat$labels == 'a1'] <-
  sample(ip.mat$values[!is.na(ip.mat$labels) & ip.mat$labels == 'a1'], 1)

#  m2.fmfit <- read.flexmirt("~/2012/sy/fm/ms-rasch-prm.txt")
# cat(deparse(round(m2.fmfit$G1$param,6)))
fmfit <- structure(c(0.941583, 1, 0, 0, 0, -0.676556, 0.758794, -0.802595,  1.28891, 0.941583, 1, 0, 0, -0.182632, 0.897435, 1.30626, NA,  NA, 0.941583, 1, 0, 0, 0, 0.177835, -1.82185, 0.005832, -0.81109,  0.941583, 1, 0, 0, 0, -1.15962, -1.229, 0.032677, 0.4922, 0.941583,  1, 0, 0.457533, 0.324595, NA, NA, NA, NA, 0.941583, 1, 0, 0,  -2.69186, -1.04012, 1.61232, NA, NA, 0.941583, 1, 0, 0, 0, -1.38231,  0.034368, -1.214, -0.648291, 0.941583, 1, 0, 0, 0, -1.85655,  -1.17135, -0.262079, -0.531158, 0.941583, 1, 0, 0, 0, -1.29475,  -0.376539, 0.02024, 0.135187, 0.941583, 1, 0, 0, 0, -1.38279,  0.524151, -0.508742, 0.633671, 0.941583, 1, 0, 0, 0, -0.979595,  -0.048528, 0.659669, 0.544857, 0.941583, 1, 0, 0, 0, -2.09039,  -1.45472, -0.472137, -0.666386, 0.941583, 1, 0, 0, 0.174682,  0.645437, 0.907132, NA, NA, 0.941583, 1, 0, 0, -0.842216, 0.490717,  1.28034, NA, NA, 0.941583, 1, 0, 0, 0, -0.913355, -0.319602,  -0.310164, -0.15536, 0.941583, 1, 0, 0, 0, 0.567085, -1.56762,  0.884553, 0.122113, 0.941583, 1, 0, 0, 0, -0.152985, -0.341317,  -0.183837, 1.17952, 0.941583, 1, 0, 0, 0, 0.168869, -0.490354,  0.373892, 1.29714, 0.941583, 1, 0, 0, 0, -0.827385, 0.626197,  -1.52994, 0.494209, 0.941583, 1, 0, 0, 0, 0.511263, -0.750358,  1.01852, 0.840026, 0.941583, 1, 0, 0, 0, 0.968905, -0.009671,  1.52297, 1.69255, 0.941583, 1, 0, 0, 0, 1.89582, 0.051828, 2.25758,  1.52469), .Dim = c(9L, 22L), .Dimnames = list(NULL, c("i1", "i2",  "i3", "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12",  "i13", "i14", "i15", "i16", "i17", "i18", "i19", "i20", "i21",  "i22")))
#  ip.mat$values <- m2.fmfit$G1$param

if (1) {
  cip.mat <- ip.mat
  cip.mat$values <- fmfit
  cM <- mxModel(model="ms", cip.mat,
                mxData(observed=m2.data, type="raw"),
                mxExpectationBA81(ItemSpec=m2.spec),
                mxFitFunctionML(),
		mxComputeOnce('fitfunction', 'fit'))
  cM <- mxRun(cM, silent=TRUE)
  omxCheckCloseEnough(cM$fitfunction$result, 50661.38, .01)
}

plan <- mxComputeSequence(steps=list(
  mxComputeEM('expectation', 'scores',
              mxComputeNewtonRaphson(freeSet='item', verbose=0L),
              information="mr1991", infoArgs=list(fitfunction='fitfunction')),
  mxComputeStandardError(),
  mxComputeHessianQuality()))

m2 <- mxModel(model="m2", ip.mat,
              mxData(observed=m2.data, type="raw"),
              mxExpectationBA81(ItemSpec=m2.spec),
              mxFitFunctionML(),
              plan)
#  m2 <- mxOption(m2, "Number of Threads", 1)
m2 <- mxRun(m2, silent=TRUE)
omxCheckCloseEnough(m2$output$minimum, 50661.377, .01)

omxCheckCloseEnough(log(m2$output$conditionNumber), 6.57, .5)
#omxCheckTrue(is.na(m2$output$conditionNumber))
#cat(deparse(round(c(m2$output$standardErrors), 3)))

semse <- c(0.022, 0.095, 0.116, 0.116, 0.108, 0.176, 0.222, 0.305, 0.382,  0.359, 0.244,
           0.215, 0.105, 0.082, 0.067, 0.07, 0.185, 0.215,  0.134, 0.061, 0.071, 0.25,
           0.244, 0.231, 0.155, 0.328, 0.209,  0.177, 0.16, 0.211, 0.176, 0.182, 0.185,
           0.187, 0.189, 0.201,  0.194, 0.174, 0.161, 0.2, 0.234, 0.409, 0.236, 0.179,
           0.154,  0.064, 0.078, 0.092, 0.084, 0.074, 0.092, 0.584, 0.493, 0.441,  0.362,
           0.1, 0.097, 0.079, 0.085, 0.113, 0.115, 0.102, 0.111,  0.079, 0.082, 0.076,
           0.092, 0.541, 0.607, 0.554, 0.337, 0.081,  0.083, 0.083, 0.098, 0.072, 0.084,
           0.103, 0.138, 0.084, 0.103,  0.141, 0.178)
#max(abs(c(m2$output$standardErrors) - semse))
omxCheckCloseEnough(c(m2$output$standardErrors), semse, .01) # similar to flexMIRT

emstat <- m2$compute$steps[[1]]$output
omxCheckCloseEnough(emstat$EMcycles, 19, 2)
omxCheckCloseEnough(emstat$totalMstep, 73, 10)
omxCheckCloseEnough(emstat$semProbeCount / length(semse), 3, .1)
omxCheckCloseEnough(m2$output$evaluations, 1062, 5)

#print(m2$matrices$item$values - fmfit)
print(m2$output$backendTime)

n <- apply(!is.na(m2.data), 2, sum)

i1 <- mxModel(m2,
              mxComputeSequence(steps=list(
                mxComputeOnce('fitfunction', 'information', "meat"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i1 <- mxRun(i1, silent=TRUE)

omxCheckTrue(i1$output$infoDefinite)
omxCheckCloseEnough(log(i1$output$conditionNumber), 7.3, .5)

#cat(deparse(round(c(i1$output$standardErrors), 3)))
se <- c(0.019, 0.1, 0.123, 0.121, 0.119, 0.237, 0.246, 0.33, 0.417,  0.386, 0.281,
        0.24, 0.108, 0.086, 0.072, 0.076, 0.221, 0.265,  0.138, 0.068, 0.085, 0.275,
        0.267, 0.263, 0.196, 0.359, 0.237,  0.208, 0.203, 0.227, 0.191, 0.199, 0.225,
        0.21, 0.215, 0.232,  0.235, 0.184, 0.179, 0.218, 0.254, 0.437, 0.26, 0.201, 0.194,
        0.07, 0.083, 0.101, 0.089, 0.079, 0.096, 0.649, 0.549, 0.507,  0.421, 0.106, 0.102,
        0.084, 0.093, 0.125, 0.124, 0.112, 0.127,  0.088, 0.089, 0.087, 0.109, 0.633, 0.704,
        0.61, 0.415, 0.089,  0.089, 0.09, 0.112, 0.083, 0.092, 0.115, 0.17, 0.095, 0.11, 0.16,  0.192)
omxCheckCloseEnough(c(i1$output$standardErrors), se, .001)  # matches flexmirt

if (0) {
  library(mirt)
  rdata <- sapply(m2.data, unclass)-1
  # for flexMIRT, write CSV
  #write.table(rdata, file="ifa-drm-mg.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  pars <- mirt(rdata, 1, itemtype="Rasch", D=1, quadpts=49, pars='values')
#  pars[pars$name=="a1",'value'] <- 1
#  pars[pars$name=="a1",'est'] <- FALSE
#  pars[pars$name=="COV_11",'est'] <- TRUE
  fit <- mirt(rdata, 1, itemtype="Rasch", D=1, quadpts=49, pars=pars, SE=TRUE, SE.type="crossprod")
  # LL -25330.691 * -2 = 50661.38
  got <- coef(fit)
}
