#options(error = browser)
require(OpenMx)
require(rpf)

mcar <- function(ret, mcar) {
  size <- prod(dim(ret))
  mask <- rep(FALSE, size)
  mask[sample.int(size, size * mcar)] <- TRUE
  shaped.mask <- array(mask, dim=dim(ret))
  ret[shaped.mask] <- NA
  ret
}

set.seed(1)

numItems <- 12
items <- list()
items[1:numItems] <- list(rpf.grm())
correct.mat <- sapply(items, rpf.rparam, version=1)

slicen <- 50
data <- rpf.sample(slicen, items, correct.mat)
for (m in seq(.1, .9, .05)) {
  data <- rbind(data, mcar(rpf.sample(slicen, items, correct.mat), m))
}

dimnames(correct.mat) <- list(c('f1', 'b'), colnames(data))

result <- expand.grid(ips=0:numItems, v=NA)

for (r in 1:nrow(result)) {
  grp <- list(spec=items,
              param=correct.mat,
              data=data,
              minItemsPerScore=result$ips[r])

  sc <- EAPscores(grp)
  
  v <- var(sc[,'f1'], na.rm=TRUE)
  result$v[r] <- v
}

c1 <- coef(lm(v ~ ips, result))
#cat(deparse(round(c1,4)))
omxCheckCloseEnough(c1, c(0.5763, 0.0144), .001)

# ------------------------------

grp <- list(spec=items,
            param=correct.mat,
            data=data,
            minItemsPerScore=ncol(correct.mat)+1)

omxCheckError(EAPscores(grp),
	      "minItemsPerScore (=13) cannot be larger than the number of items (=12)")
