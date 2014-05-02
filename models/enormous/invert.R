library(OpenMx)
library(rbenchmark)
library(ggplot2)
library(reshape2)

# create random symmetric positive definite matrix
genSymm <- function (n, ev = runif(n, 1, 10)) {
  if (n==1) return(matrix(ev^2,1,1))
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  if (0) {
    toneg <- is.na(match(1:n, sample.int(n, n/2)))
    ev[toneg] <- -ev[toneg]
  }
  Z <- t(O) %*% diag(ev) %*% O
  Z
}

nonZero <- function(mat) sum(mat!=0) / prod(dim(mat))

mkSparse <- function(n, sparse) {
  mat <- diag(exp(runif(n, -3, 6)))
  while(nonZero(mat) < sparse) {
    bsize <- sample.int(n/2, 1)
    block <- genSymm(bsize)
    corner <- sample.int(n - bsize, 1)
    bottom <- corner + bsize - 1
    mat[corner:bottom, corner:bottom] <- mat[corner:bottom, corner:bottom] + block
  }
  mat
}

grid <- expand.grid(size=seq(10,1000,length.out=5),
                    sparse=seq(.01, .2, length.out=5),
                    dtm=NA, stm=NA)

schulz1933 <- function(orig) {
  iter <- 0
  I <- diag(nrow(orig))
  prev <- I / norm(orig, 'f')
  AV <- orig %*% prev
  while (iter < 100) {
    iter <- iter + 1
    nguess <- prev %*% (2 * I - AV)
    AV <- orig %*% nguess
    err <- norm(AV-I, "1")
    if (err < 1e-6) return(list(iter=iter, out=nguess))
    if (err > 1) stop(paste("schulz1933 failed", err))
    prev <- nguess
  }
}

soleymani2013 <- function(orig) {
  iter <- 0
  I <- diag(nrow(orig))
  prev <- I / norm(orig, 'f')
  AV <- orig %*% prev
  while (iter < 100) {
    iter <- iter + 1
    rpart <- 7*I+AV %*% (-21*I + AV %*% (35 * I + AV %*% (-35 * I + AV %*% (21 * I + AV %*% (-7*I + AV)))))
    nguess <- prev %*% rpart
    AV <- orig %*% nguess
    err <- norm(AV-I, "2")
    if (err < 1e-6) return(list(iter=iter, out=nguess))
    if (err > 1) stop(paste("soleymani2013 failed", err))
    prev <- nguess
  }
}

set.seed(1)
for (gx in 1:nrow(grid)) {
  print(c(grid$size[gx], grid$sparse[gx]))
  mat <- mkSparse(grid$size[gx], grid$sparse[gx])
#  schulz1933(mat)
  df <- benchmark(dense=solve(mat),
#                  sparse=soleymani2013(mat),
                  sparse=imxSparseInvert(mat),
                  replications=2)
  rownames(df) <- df$test
  grid$dtm[gx] <- df['dense','elapsed']
  grid$stm[gx] <- df['sparse','elapsed']
  print(gx)
}

pdat <- melt(grid, id.vars=c("size","sparse"), variable.name="algorithm", value.name="time")

ggplot(pdat, aes(size, time, color=algorithm)) + geom_line() + facet_wrap(~sparse)
