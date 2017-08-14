# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
library(mvtnorm)

set.seed(1)

numItems <- 16
numPeople <- 1000

items <- list()
correct <- list()
items[1:numItems] <- list(rpf.grm(factor=6))
correct <- sapply(items, rpf.rparam, version=1)
correct['a1', 12:16] <- 0
correct['a2', 1:5] <- 0
correct['a3', 5:16] <- 0
correct['a4', c(1:4, 9:16)] <- 0
correct['a5', c(1:8, 13:16)] <- 0
correct['a6', 1:12] <- 0

maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))

true.mean <- c(.2, -.1, .3, -.2, -.5, .15)
true.cov <- diag(6)
diag(true.cov) <- c(.5, 1.2, 1.5, .9, 1.1, 1.3)
true.cov[1,2] <- -.2
true.cov[2,1] <- -.2
true.pt <- c(true.mean, diag(true.cov)[-2], true.cov[2,1])

data <- rpf.sample(numPeople, items, correct, mean=true.mean, cov=true.cov)

ip.mat <- mxMatrix(name="item", nrow=maxParam, ncol=numItems,
                   values=correct, free=FALSE)
rownames(ip.mat) <- c(paste('f', 1:6, sep=""), rep('p', nrow(ip.mat)-6))
colnames(ip.mat) <- colnames(data)
col6mask <- correct[,6] != 0
ip.mat$free[col6mask,6] <- TRUE

cache <- list(c(-188.63842, -303.37294, 562.09478, 715.55771, 405.68351, 732.46805,  
                -160.01455, 1884.34914, 98.36976, 73.19764, 124.96415, -329.53384,  
                -29.48787, -1331.09128, -14.59534, 15.53051, -80.80697),
              c(-93.28779, 188.45762, -233.45221, 384.96591, -186.08688, 581.75515,  
                -22.35932, -157.65153, 465.58903, -88.02069, 143.7526, -81.11557,  -13.56995, -50.40794, -56.54828, 2.67365, 244.76936),
              c(109.12714, 212.34538, -195.441, 391.51011, 1241.95121, 2148.99911,  
                -260.00634, 729.12616, 877.78973, 145.0623, -328.81856, -1977.33925,  
                -1.32911, -158.44536, -165.04632, -41.28467, -1927.31541),
              c(31.12355, 52.74843, -115.22288, -183.10723, -64.0352, 107.29031,  
                51.47168, -255.91144, 734.44948, 148.68402, 171.971, 73.33686,  7.71757, -0.75168, -119.34241, 4.36869, 67.36673))

trials <- 4
diff <- matrix(NA, trials, 17)

for (seed in 1:trials) {
  set.seed(seed)
  
  ip.mat$values[col6mask,6] <- rpf.rparam(items[[6]], version=1)[col6mask]
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=6, values=runif(6, -1, 1), free=TRUE)
  colnames(m.mat) <- paste('f', 1:6, sep="")
  var <- runif(2, .1, 3)
  cov1 <- runif(1, -min(abs(var)), min(abs(var)))
  cov <- diag(runif(6, .1, 3))
  cov[1:2,1:2] <- matrix(c(var[1],cov1,cov1,var[2]), nrow=2)
  
  cov.mat <- mxMatrix(name="cov", nrow=6, ncol=6, values=cov)
  dimnames(cov.mat) <- list(paste('f', 1:6, sep=""), paste('f', 1:6, sep=""))
  cov.mat$free <- cov.mat$values != 0
  cov.mat$labels[1:2,1:2] <- matrix(c("v1","c12","c12","v2"), nrow=2)

  m1 <- mxModel(model="latent",
                ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=items, mean="mean", cov="cov",
                                  qwidth=5, qpoints=21),
                mxFitFunctionML(),
                mxComputeSequence(steps=list(
                  mxComputeOnce('fitfunction', 'gradient'),
                  mxComputeReportDeriv())))
  
  objective1 <- function(x) {
    ip.mat$values[col6mask,6] <- x[1:4]
    m.mat$values[1:6] <- x[5:10]
    diag(cov.mat$values) <- x[11:16]
    cov.mat$values[1,2] <- x[17]
    cov.mat$values[2,1] <- x[17]
    m1.probe <- mxModel(m1, m.mat, cov.mat, ip.mat,
                        mxComputeSequence(steps=list(
                          mxComputeOnce('fitfunction', 'fit'))))
    m1.probe <- mxRun(m1.probe, silent=TRUE)
    got <- m1.probe$output$fit
    got
  }
  
  probe.pt <- c(ip.mat$values[col6mask,6], m.mat$values, diag(cov.mat$values), cov.mat$values[2,1])
#  print(probe.pt - true.pt)
  
  deriv <- c()
  if (0) {
    require(numDeriv)
    deriv <- grad(objective1, probe.pt, method="simple")
    cat(deparse(round(deriv,5)), fill=TRUE)
  } else {
    deriv <- cache[[seed]]
  }

  m1 <- mxRun(m1, silent=TRUE)
  #  print(probe.pt)
#print(m1$output$gradient)
  covTerm <- match('c12', names(m1$output$gradient))
  aGrad <- c(m1$output$gradient[-covTerm], m1$output$gradient[covTerm])
  
  diff[seed,] <- aGrad - deriv
}
#round(apply(abs(diff), 2, median),2)
print(max(apply(abs(diff), 2, median)))
omxCheckCloseEnough(apply(abs(diff), 2, median), rep(0,dim(diff)[2]), 6.1)
