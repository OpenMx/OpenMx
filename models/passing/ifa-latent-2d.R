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
	correct[[ix]] <- rpf.rparam(items[[ix]])
}
correct.mat <- simplify2array(correct)

maxParam <- max(vapply(items, function(i) rpf.numParam(i), 0))

true.mean <- c(.2,-.1)
true.cov <- matrix(c(.5, -.2, -.2, .5), nrow=2)
data <- rpf.sample(numPeople, items, correct.mat, mean=true.mean, cov=true.cov)

ip.mat <- mxMatrix(name="ItemParam", nrow=maxParam, ncol=numItems,
                   values=correct.mat, free=FALSE)

cache <- list(c(-386.18928, -252.51883, 35.11062, 74.30228, 55.13235),
              c(-218.59169, 146.25932, 61.66015, 243.46282, 2.52876),
              c(-251.25287, 143.21162, 61.76765, 251.09119, 43.17802),
              c(-117.41808, -380.62556, 69.23508, 197.80985, -14.66206),
              c(-152.39917, 55.13713, 96.17968, 160.79444, 44.92327),
              c(155.65986, 365.33715, 57.25361, 173.34409, -14.9723),
              c(1150.23715, 809.20463, -852.3527, -1269.55056, -506.94031),
              c(-188.0905, -237.78866, 86.86674, 116.09097, 50.44799),
              c(-993.15851, -927.75109, -405.12654, -709.25762, -349.15984),
              c(-299.16442, -267.18778, 41.78237, 53.74096, 41.61348))

trials <- 10
diff <- matrix(NA, trials, 5)
for (seed in 1:trials) {
  set.seed(seed)
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=2, values=runif(2, -1, 1), free=TRUE)
  var <- runif(2, .1, 3)
  cov1 <- runif(1, -min(abs(var)), min(abs(var)))
  cov <- matrix(c(var[1],cov1,cov1,var[2]), nrow=2)
  
  cov.mat <- mxMatrix(name="cov", nrow=2, ncol=2, values=cov,
                      free=TRUE, labels=c("v1","c12","c12","v2"))

  m1 <- mxModel(model="latent",
                ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=items,
                                  ItemParam="ItemParam"),
                mxFitFunctionML(),
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation'),
                  mxComputeOnce('fitfunction', gradient=TRUE))))
  
  objective1 <- function(x) {
    m.mat@values[1:2] <- x[1:2]
    cov.mat@values[1,1] <- x[3]
    cov.mat@values[1,2] <- x[4]
    cov.mat@values[2,1] <- x[4]
    cov.mat@values[2,2] <- x[5]
    m1.probe <- mxModel(m1, m.mat, cov.mat,
                        mxComputeSequence(steps=list(
                          mxComputeOnce('expectation'),
                          mxComputeOnce('fitfunction', fit=TRUE))))
    m1.probe <- mxRun(m1.probe, silent=TRUE)
    got <- m1.probe@output$minimum
    #  print(got)
    got
  }
  
  probe.pt <- c(m.mat@values, cov.mat@values[1:2,1], cov.mat@values[2,2])

  deriv <- c()
  if (0) {
    require(numDeriv)
    deriv <- grad(objective1, probe.pt, method="simple")
    cat(deparse(round(deriv,5)), fill=TRUE)
  } else {
    deriv <- cache[[seed]]
  }
  
  m1 <- mxRun(m1, silent=TRUE)
  diff[seed,] <- m1@output$gradient - deriv
#  print(probe.pt)
}
#apply(abs(diff), 2, median)
omxCheckCloseEnough(apply(abs(diff), 2, median), rep(0,5), .06)
