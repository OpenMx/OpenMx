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
ip.mat$free[,6] <- TRUE

cache <- list(c(77.55278, -73.91341, -166.20929, 182.00042, -0.95405, -481.91857,  
                -411.68003, -224.00015, 579.44648, -502.315, 84.53358, -42.47443,  -24.27724, -32.7925, -155.6535, -80.4796, 148.63529),
              c(-30.00643, -42.32397, -30.12896, -322.33044, -254.95183, 137.4676,  
                -382.20775, -542.30859, 785.89643, -277.49958, -79.97633, -1.05163,  
                -356.0611, -83.50646, -171.87116, -29.25338, 361.60296),
              c(-178.51635, -1.93742, -120.78015, -481.66085, 228.42419, 435.21254,  
                -169.77725, -274.75197, 867.96436, -105.45071, -38.42063, -113.69795,  
                -72.34722, -138.71124, -539.72992, -100.47599, -312.49079),
              c(15.86328, -10.24639, 8.74378, 17.08915, 730.40302, -759.16442,  
                436.34736, -209.00355, -242.96903, -134.40352, -238.60034, -273.59263,  
                -71.49758, -21.53199, 7.95964, -39.02151, 789.94771),
              c(-222.46406, 162.36542, -109.31309, -242.00755, 272.85096, -1019.65358,  
                -186.77005, -18.78546, 21.45854, -484.56252, 15.01884, -346.89958,  22.43994, 1.49772, 25.12492, -152.2727, 115.26171),
              c(-5.20869, -130.04494, -43.9493, 255.75791, -17.93615, -1302.40378,  
                -78.27092, 352.63914, -297.34807, -936.77417, 94.79776, -597.26737,  7.17598, -56.73478, 1.91608, -356.48808, 74.11117),
              c(40.13739, -19.16884, -73.44294, -176.43786, -1120.72168, -85.27712,  
                -1003.15307, -270.76263, 891.44778, -379.24816, -116.80423, -231.10691,  
                -342.16178, -25.94346, -266.95926, -51.34715, -73.11103),
              c(-180.84729, 41.56373, 137.62677, -266.77822, 215.6258, 117.36318,  
                -72.02323, -260.49811, 489.90572, 91.92744, -58.75196, 108.71605,  0.55119, -74.87482, -70.4871, -20.43705, -18.71787))

trials <- 8
diff <- matrix(NA, trials, 17)

for (seed in 1:trials) {
  set.seed(seed)
  
  ip.mat$values[,6] <- rpf.rparam(items[[6]])
  
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
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=items, design=design,
                                  ItemParam="ItemParam", qwidth=5, qpoints=21),
                mxFitFunctionML(),
                mxComputeSequence(steps=list(
                  mxComputeOnce('fitfunction', 'gradient'),
                  mxComputeReportDeriv())))
  
  objective1 <- function(x) {
    ip.mat$values[,6] <- x[1:4]
    m.mat$values[1:6] <- x[5:10]
    diag(cov.mat$values) <- x[11:16]
    cov.mat$values[1,2] <- x[17]
    cov.mat$values[2,1] <- x[17]
    m1.probe <- mxModel(m1, m.mat, cov.mat, ip.mat,
                        mxComputeSequence(steps=list(
                          mxComputeOnce('fitfunction', 'fit'))))
    m1.probe <- mxRun(m1.probe, silent=TRUE)
    got <- m1.probe$output$minimum
    #  print(got)
    got
  }
  
  probe.pt <- c(ip.mat$values[,6], m.mat$values, diag(cov.mat$values), cov.mat$values[2,1])
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
  covTerm <- match('c12', names(m1$output$gradient))
  aGrad <- c(m1$output$gradient[-covTerm], m1$output$gradient[covTerm])
  
  diff[seed,] <- aGrad - deriv
}
#round(apply(abs(diff), 2, median),2)
omxCheckCloseEnough(apply(abs(diff), 2, median), rep(0,dim(diff)[2]), 3)
