library(mvtnorm)
library(OpenMx)

set.seed(1)

cov <- matrix(0, 12, 12)
cov[1:4,1:4] <- rWishart(1, 4, diag(4))[,,1]
cov[5:8,5:8] <- rWishart(1, 4, diag(4))[,,1]
cov[9:12,9:12] <- rWishart(1, 4, diag(4))[,,1]

mean <- rnorm(12, sd=sqrt(diag(cov)))

mxOption(NULL, "maxOrdinalPerBlock", 12)
lk1 <- omxMnor(cov, mean, matrix(-1, 12, 1), matrix(1, 12, 1))
omxCheckCloseEnough(lk1, 1.41528651675062e-05, 1e-7)

mxOption(NULL, "maxOrdinalPerBlock", 4)
lk2 <- omxMnor(cov, mean, matrix(-1, 12, 1), matrix(1, 12, 1))
omxCheckCloseEnough(lk1, lk2, 1e-7)

mxOption(NULL, "maxOrdinalPerBlock", 3)
lk3 <- omxMnor(cov, mean, matrix(-1, 12, 1), matrix(1, 12, 1))
omxCheckTrue(lk1 != lk3)
omxCheckCloseEnough(lk1, lk3, 5e-6)

mxOption(NULL, "maxOrdinalPerBlock", 2)
lk4 <- omxMnor(cov, mean, matrix(-1, 12, 1), matrix(1, 12, 1))
omxCheckTrue(lk1 != lk4)
omxCheckCloseEnough(lk1, lk4, 1e-5)

mxOption(NULL, "maxOrdinalPerBlock", 1)
lk5 <- omxMnor(cov, mean, matrix(-1, 12, 1), matrix(1, 12, 1))
omxCheckTrue(lk1 != lk5)
omxCheckCloseEnough(lk1, lk5, 1e-4)

# ----------------

cov <- diag(rlnorm(2))
mean <- matrix(runif(2), 2, 1)

mxOption(NULL, "maxOrdinalPerBlock", 2)
lk1 <- omxMnor(cov, mean, matrix(c(-1,-Inf), 2, 1), matrix(c(Inf,1), 2, 1))
omxCheckCloseEnough(lk1,
                    pmvnorm(lower=c(-1,-Inf), upper=c(Inf,1),
                            mean=c(mean), sigma=cov))

mxOption(NULL, "maxOrdinalPerBlock", 1)
lk2 <- omxMnor(cov, mean, matrix(c(-1,-Inf), 2, 1),
               matrix(c(Inf,1), 2, 1))
omxCheckCloseEnough(lk1, lk2)

omxCheckEquals(omxMnor(cov, mean,
                       matrix(c(-Inf,-Inf), 2, 1),
                       matrix(c(Inf,Inf), 2, 1)), 1.0)

# ----------------

blocks <- 10
perBlock <- 5

cov <- matrix(0, blocks*perBlock, blocks*perBlock)
for (bl in 1:blocks) {
  ind <- seq(1+(bl-1)*perBlock, bl*perBlock)
  cov[ind, ind] <- rWishart(1, perBlock*2, diag(perBlock))[,,1]
}

mean <- rnorm(nrow(cov), sd=sqrt(diag(cov)))

mxOption(NULL, "maxOrdinalPerBlock", 12)
lk1 <- omxMnor(cov, mean,
               matrix(-1, blocks*perBlock, 1),
               matrix(1, blocks*perBlock, 1))

omxCheckCloseEnough(log(lk1), -115.15, .1)

