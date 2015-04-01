# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
library(mvtnorm)

set.seed(1)
numItems <- 12
numPeople <- 500

items <- list()
correct <- list()
for (ix in 1:numItems) {
	items[[ix]] <- rpf.grm(factors=2)
	correct[[ix]] <- rpf.rparam(items[[ix]], version=1)
}
correct.mat <- simplify2array(correct)

maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))

true.mean <- c(.2,-.1)
true.cov <- matrix(c(.5, -.2, -.2, .5), nrow=2)
data <- rpf.sample(numPeople, items, correct.mat, mean=true.mean, cov=true.cov)

ip.mat <- mxMatrix(name="item", nrow=maxParam, ncol=numItems,
                   values=correct.mat, free=FALSE)
rownames(ip.mat) <- c(paste('f', 1:2, sep=""), rep('p', nrow(ip.mat)-2))
colnames(ip.mat) <- colnames(data)
ip.mat$free[,1:2] <- TRUE

cache <- list(c(32.64709, 40.19653, 138.53231, -39.95988, 55.5181, 70.53231,  -327.55366, -217.73847, 55.02332, 107.86482, 65.56429),
              c(7.09441, 89.77065, 154.90115, -27.91789, 49.86673, 86.90115,  -166.89402, 146.50368, 75.37147, 219.22294, 22.90383),
              c(-16.86464, 112.33415, 156.67213, -45.87323, 69.06868, 88.67213,  -178.78624, 154.45145, 82.23031, 233.08362, 50.29965),
              c(116.77814, -49.40159, 144.25094, 37.88699, -30.23027, 76.25094,  
                -110.05324, -332.33857, 76.33498, 186.92474, 13.85316),
              c(57.40723, 24.16563, 155.26433, -23.39773, 35.98561, 87.26433,  -106.70071, 92.57732, 102.43476, 169.8306, 50.8937),
              c(15.29993, 102.66603, 182.04623, -19.57674, 56.71805, 114.04623,  210.79779, 382.99158, 56.71743, 128.95432, -12.30778),
              c(163.98019, -74.3645, 242.23824, 79.13552, -52.32248, 174.23824,  
                1216.51275, 959.54392, -942.59298, -1528.70604, -637.4546),
              c(100.77217, -21.86135, 145.69112, -3.81738, 20.32068, 77.69112,  -160.5544, -202.25654, 89.84531, 134.02613, 60.17379))

trials <- 8
diff <- matrix(NA, trials, 11)
for (seed in 1:trials) {
  set.seed(seed)
  
  ip.mat$values[,1:2] <- c(1.1,1,0)
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=runif(2, -1, 1), free=TRUE)
  colnames(m.mat) <- paste('f', 1:2, sep="")
  var <- runif(2, .1, 3)
  cov1 <- runif(1, -min(abs(var)), min(abs(var)))
  cov <- matrix(c(var[1],cov1,cov1,var[2]), nrow=2)
  
  cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=cov,
                      free=TRUE, labels=c("v1","c12","c12","v2"))
  dimnames(cov.mat) <- list(paste('f', 1:2, sep=""), paste('f', 1:2, sep=""))

  m1 <- mxModel(model="latent",
                ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=items, mean="mean", cov="cov"),
                mxFitFunctionML(),
                mxComputeSequence(list(
		    mxComputeOnce('fitfunction', 'gradient'),
        mxComputeReportDeriv())))
  
  objective1 <- function(x) {
    ip.mat$values[,1:2] <- x[1:6]
    m.mat$values[1:2] <- x[7:8]
    cov.mat$values[1,1] <- x[9]
    cov.mat$values[1,2] <- x[10]
    cov.mat$values[2,1] <- x[10]
    cov.mat$values[2,2] <- x[11]
    m1.probe <- mxModel(m1, ip.mat, m.mat, cov.mat,
                        mxComputeSequence(steps=list(
                          mxComputeOnce('fitfunction', 'fit'))))
    m1.probe <- mxRun(m1.probe, silent=TRUE)
    got <- m1.probe$output$minimum
    #  print(got)
    got
  }
  
  probe.pt <- c(ip.mat[,1:2], m.mat$values, cov.mat$values[1:2,1], cov.mat$values[2,2])

  deriv <- c()
  if (0) {
    require(numDeriv)
    deriv <- grad(objective1, probe.pt, method="simple")
    cat(deparse(round(deriv,5)), fill=TRUE)
  } else {
    deriv <- cache[[seed]]
  }
  
  m1 <- mxRun(m1, silent=TRUE)
  diff[seed,] <- m1$output$gradient - deriv
#  print(probe.pt)
}
#apply(abs(diff), 2, median)
omxCheckCloseEnough(apply(abs(diff), 2, median), rep(0,11), .1)
