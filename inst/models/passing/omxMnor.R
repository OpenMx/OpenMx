#
#   Copyright 2007-2018 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


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
foo <- try(omxMnor(cov, mean, matrix(-1, 12, 1), matrix(1, 12, 1)))
omxCheckTrue(is(foo, 'try-error'))
omxCheckTrue(foo[1] == "Error in omxMnor(cov, mean, matrix(-1, 12, 1), matrix(1, 12, 1)) : \n  Ordinal covariance has dependent block larger than 3x3. You must increase mxOption maxOrdinalPerBlock\n")

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

