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
  mxPenaltyLASSO(what = paste0('c', 1:(D * (D-1)/2))[regGroup=='lasso'],
                    name = "lasso", epsilon=.01),
  mxPenaltyRidge(what = paste0('c', 1:(D * (D-1)/2))[regGroup=='ridge'],
                         name = "ridge", epsilon=.01),
  mxPenaltyElasticNet(what = paste0('c', 1:(D * (D-1)/2))[regGroup=='elasticNet'],
                    name = "en", epsilon=.01, lambda.step=.01,
                    alpha = .4, alpha.max = .4))

m1 <- mxModel(m1, lapply(
  m1$penalties,
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

expect_equal(mxEval(lasso, m1), matrix(0,0,0))
expect_equivalent(mxEval(lasso, m1, compute = TRUE), matrix(0.4, 1, 1))

for (rep in 1:10) {
  fit1 <- mxModel(m1,
                  mxData(randomCov(), type = 'cov', numObs = 100))


  mxOption(key="Analytic Gradients", value="Yes")

  fit1 <- mxRun(mxModel(fit1,
                        mxComputeSequence(list(
                          mxComputeOnce('fitfunction', c('fit','gradient')),
                          mxComputeReportDeriv()))), silent = TRUE)
  expect_equal(fit1$output$evaluations, 2)

  # Added correctly?  
  expect_equal(sum(
    fit1$fitfunction$result[1,1],
    sapply(fit1$penalties, function(x) x$result[1,1])),
    fit1$output$fit)

  # Also test mxEval
  expect_equal(fit1$fitfunction$result[1,1],
               mxEval(fitfunction, fit1, compute = TRUE)[1,1])
  for (p1 in fit1$penalties) {
    expect_equivalent(fit1[[p1$name]]$result[1,1],
                      mxEvalByName(p1$name, fit1, compute = TRUE)[1,1])
  }

  mxOption(key="Analytic Gradients", value="No")

  fit2 <- mxRun(mxModel(fit1,
                        mxComputeSequence(list(
                          mxComputeOnce('fitfunction', c('fit','gradient')),
                          mxComputeReportDeriv()))), silent = TRUE)
  expect_true(fit2$output$evaluations > 2)

  expect_equal(fit1$output$gradient, fit2$output$gradient, 1e-6)
}

m1 <- expect_warning(mxPenaltySearch(m1),
                     "model does not satisfy the first-order optimality conditions")
expect_equal((-2*logLik(m1))[1],m1$output$fit)
if(mxOption(NULL,"Default optimizer") != "CSOLNP"){
	expect_equal((-2*logLik(m1))[1],876.4839,5e-5)
}

detail <- m1$compute$steps$PS$output$detail

limit <- c(lasso=1, ridge=11, elasticNet=2)

for (cx in 1:(D * (D-1)/2)) {
  col <- paste0('c', cx)
  #print(paste(regGroup[cx], col))
  val <- abs(detail[[col]])
  val[val < 0.01] <- 0
  if (all(diff(val) <= 0)) next
  thr <- limit[regGroup[cx]]
  expect_equivalent(table(diff(val) <= 0)["FALSE"], 0, thr)
}

# Try with fixed hyperparam
m3 <- omxSetParameters(m1, labels = c('lambda','alpha'), free = FALSE)
expect_equivalent(mxEval(lasso, m3),
                  mxEval(lasso, m3, compute = TRUE))
m3 <- mxRun(m3)
expect_equivalent(mxEval(lasso, m3),
                  mxEval(lasso, m3, compute = TRUE))

m2 <- mxPenaltyZap(m1, silent = TRUE)
expect_error(summary(m2), "This model is in an inconsistent state.")
fit2 <- mxRun(m2)
expect_equal(fit2$output$fit, 853.759, 1e-2)

m3 <- m1
m3$fitfunction$applyPenalty <- FALSE
fit3 <- mxRun(m3)

expect_equivalent(sapply(fit3$penalties, function(x) x$result[1,1]),
             rep(0,3))
expect_equal(fit3$fitfunction$result[1,1], fit3$output$fit)
expect_equal(fit3$output$fit, 822.757, 5e-5)
expect_true(mxEval(lasso, fit3, compute = TRUE) != 0)

#Sanity check:
if(mxOption(NULL,"Default optimizer") != "CSOLNP"){
	expect_equal(
		as.vector(coef(m1)),
		c(1.24,-0.25,1.72,-0.01,0.74,1.72,0.01,0.69,-0.37,1.86,0.61,0.09,-0.26,0.53,1.62,-0.44,0.29,-0.32,0.01,0.26,1.89,12.0,0.40),
		1e-2)
}
