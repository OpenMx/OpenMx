library(OpenMx)
library(testthat)

set.seed(1)

D <- 6

randomCov <- function() {
  cf <- matrix(0, D,D)
  cf[lower.tri(cf, FALSE)] <- runif(D * (D-1)/2, -.8, .8)
  diag(cf) <- rnorm(D, mean = 1.2, sd = .1)
  dat <- cf %*% t(cf)
  dimnames(dat) <- list(paste0('x',1:D),paste0('x',1:D))
  dat
}

label <- matrix(NA, D,D)
diag(label) <- paste0('d',1:D)
label[lower.tri(label)] <- paste0('c', 1:(D * (D-1)/2))

regGroup <- cut(1:(D * (D-1)/2), 3)
levels(regGroup) <- c('lasso', 'ridge', 'elasticNet')

m1 <- mxModel(
  'm1',
  mxData(randomCov(), type = 'cov', numObs = 100),
  mxMatrix('Symm', D,D, TRUE, values=diag(D) + .2,
           labels=label[lower.tri(label,TRUE)],
           dimnames=list(paste0('x',1:D),paste0('x',1:D)),
           name="cov"),
  mxMatrix(nrow=1, ncol=2, free=TRUE, values=c(.4, .4),
           labels=c('lambda','alpha')),
  mxExpectationNormal('cov'),
  mxFitFunctionML(),
  mxRegularizeLASSO(what = paste0('c', 1:(D * (D-1)/2))[regGroup=='lasso'],
                    name = "lasso", epsilon=.01),
  mxRegularizeRidge(what = paste0('c', 1:(D * (D-1)/2))[regGroup=='ridge'],
                         name = "ridge", epsilon=.01),
  mxRegularizeElasticNet(what = paste0('c', 1:(D * (D-1)/2))[regGroup=='elasticNet'],
                    name = "en", epsilon=.01, lambda.step=.01,
                    alpha = .4, alpha.max = .4))

m1 <- mxModel(m1, lapply(
  m1$regularizations,
  function(reg) {
    reg$hpranges$lambda <- 400*reg$hpranges$lambda
    reg
  }))

if (0) {
  m1 <- mxOption(m1,"Always Checkpoint","Yes")
  m1 <- mxOption(m1,"Checkpoint Units","evaluations")
  m1 <- mxOption(m1,"Checkpoint Count",1)
  m1 <- mxOption(m1, "Checkpoint Fullpath", "/dev/fd/2")
}

for (rep in 1:10) {
  fit1 <- mxModel(m1,
                  mxData(randomCov(), type = 'cov', numObs = 100))


  mxOption(key="Analytic Gradients", value="Yes")

  fit1 <- mxRun(mxModel(fit1,
                        mxComputeSequence(list(
                          mxComputeOnce('fitfunction', c('fit','gradient')),
                          mxComputeReportDeriv()))), silent = TRUE)
  expect_equal(fit1$output$evaluations, 2)

  mxOption(key="Analytic Gradients", value="No")

  fit2 <- mxRun(mxModel(fit1,
                        mxComputeSequence(list(
                          mxComputeOnce('fitfunction', c('fit','gradient')),
                          mxComputeReportDeriv()))), silent = TRUE)
  expect_true(fit2$output$evaluations > 2)

  expect_equal(fit1$output$gradient, fit2$output$gradient, 1e-6)
}

m1 <- mxRun(m1)

detail <- m1$compute$steps$REG$output$detail

limit <- c(lasso=1, ridge=11, elasticNet=2)
for (cx in 1:(D * (D-1)/2)) {
  col <- paste0('c', cx)
#  print(paste(regGroup[cx], col))
  val <- abs(detail[[col]])
  val[val < 0.01] <- 0
  if (all(diff(val) <= 0)) next
  thr <- limit[regGroup[cx]]
  expect_equivalent(table(diff(val) <= 0)["FALSE"], 0, thr)
}
