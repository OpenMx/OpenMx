# echo 0 > /proc/self/coredump_filter  # normally 023
# R --vanilla --no-save -f models/failing/bock-aitkin-1981.R
# R -d gdb --vanilla --no-save -f models/failing/bock-aitkin-1981.R

#options(error = browser)
require(OpenMx)
require(rpf)
library(mvtnorm)

set.seed(1)

numItems <- 16
numPeople <- 500

items <- list()
correct <- list()
for (ix in 1:numItems) {
  if (ix <= 5 || ix > numItems-5) {
    items[[ix]] <- rpf.grm(factors=2)
  } else {
    items[[ix]] <- rpf.grm(factors=3)
  }
	correct[[ix]] <- rpf.rparam(items[[ix]])
}
maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))
correct.mat <- sapply(correct, function(p) c(p, rep(NA, maxParam-length(p))))

true.mean <- c(.2, -.1, .3, -.2, -.5, .15)
true.cov <- diag(6)
diag(true.cov) <- c(.5, 1.2, 1.5, .9, 1.1, 1.3)
true.cov[1,2] <- -.2
true.cov[2,1] <- -.2
true.pt <- c(true.mean, diag(true.cov)[-2], true.cov[2,1])

design <- matrix(c(rep(1L,numItems-5), rep(NA,5),
                   rep(NA,5), rep(2L, numItems-5),
                   as.integer(kronecker(3:6,rep(1,4)))), byrow=TRUE, ncol=numItems)

data <- rpf.sample(numPeople, items, correct.mat, design=design, mean=true.mean, cov=true.cov)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=correct.mat, free=FALSE)

cache <- list(c(-154.47824, 68.38876, -141.9051, 168.81419, -163.4029, 70.22361,  
                34.2317, 30.47046, -67.53613, -20.37652, -17.7683, 3.54871, 64.94546 ),
              c(-513.74756, 265.22542, -300.76141, -252.1387, 515.74301, 134.22559,  
                -163.65775, 14.80394, -106.28184, -28.84007, -184.27976, -17.49414,  180.68775),
              c(-546.05265, 330.52021, -382.72443, -58.10493, 231.0593, 125.1684,  
                -153.4536, -40.19255, -83.96658, 10.302, -10.73235, -0.7226,  246.54888),
              c(107.54376, -315.77402, -174.51499, -71.11526, 238.98663, -303.64933,  
                -20.35021, -80.25071, -69.56831, -64.78418, -12.00931, -81.26354,  208.5762),
              c(-323.44541, 225.38161, 52.08564, -109.818, -92.97802, 263.88848,  
                -57.49684, -26.74753, -36.9511, -26.04965, 9.54675, -65.79383,  229.31144),
              c(7.90011, 442.91977, -188.27667, -51.21122, 440.24953, 147.90806,  
                55.09937, -52.71043, -8.95109, -51.11931, -117.02963, -11.38825,  36.38522),
              c(117.31755, 40.07052, -276.33104, -65.08144, 59.46788, 60.81853,  
                19.97912, 47.13061, -71.74223, 5.18766, -76.9569, -9.07506, 3.9718 ),
              c(53.18839, -97.95571, 249.06888, 52.84158, -61.07785, -30.62449,  
                0.18276, 37.84789, -220.10131, -20.13843, 7.1971, -25.77725,  105.98583),
              c(-554.56612, -440.73448, -823.2569, -86.44144, 120.04677, -317.86967,  
                -140.17828, -116.4752, -556.48678, 4.49604, -22.66227, -66.98681,  -81.5922),
              c(-15.57648, -261.64554, -168.64094, 307.14066, -150.18689, -289.77369,  
                44.92952, -47.7352, -14.03064, -111.88884, -0.19232, -70.95518,  91.73016))

trials <- 10
diff <- matrix(NA, trials, 13)

for (seed in 1:trials) {
  set.seed(seed)
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=6, values=runif(6, -1, 1), free=TRUE)
  var <- runif(2, .1, 3)
  cov1 <- runif(1, -min(abs(var)), min(abs(var)))
  cov <- diag(runif(6, .1, 3))
  cov[1:2,1:2] <- matrix(c(var[1],cov1,cov1,var[2]), nrow=2)
  
  cov.mat <- mxMatrix(name="cov", nrow=6, ncol=6, values=cov)
  cov.mat@free <- cov.mat@values != 0
  cov.mat@labels[1:2,1:2] <- matrix(c("v1","c12","c12","v2"), nrow=2)

  m1 <- mxModel(model="latent",
                ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=items, design=design,
                                  ItemParam="ItemParam", qwidth=5, qpoints=21),
                mxFitFunctionML(),
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation'),
                  mxComputeOnce('fitfunction', gradient=TRUE))))
  
  objective1 <- function(x) {
    m.mat@values[1:6] <- x[1:6]
    diag(cov.mat@values) <- x[7:12]
    cov.mat@values[1,2] <- x[13]
    cov.mat@values[2,1] <- x[13]
    m1.probe <- mxModel(m1, m.mat, cov.mat,
                        mxComputeSequence(steps=list(
                          mxComputeOnce('expectation'),
                          mxComputeOnce('fitfunction', fit=TRUE))))
    m1.probe <- mxRun(m1.probe, silent=TRUE)
    got <- m1.probe@output$minimum
    #  print(got)
    got
  }
  
  probe.pt <- c(m.mat@values, diag(cov.mat@values), cov.mat@values[2,1])
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
  aGrad <- c(m1@output$gradient[-8], m1@output$gradient[8])
  
  diff[seed,] <- aGrad - deriv
}
#apply(abs(diff), 2, median)
omxCheckCloseEnough(apply(abs(diff), 2, median), rep(0,13), 2.6)
