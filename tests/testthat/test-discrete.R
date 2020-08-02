library(testthat)
context("discrete")
library(OpenMx)

# to test:
# check equivalence of parameterizations (total var = 1 vs loading = 1) ??

verifyFrontBackMatch <- function(m1) {
  t1 <- mxGetExpected(m1, "thresholds")
  
  m1 <- mxRun(mxModel(m1, mxComputeSequence(list(
    mxComputeOnce('expectation'),
    mxComputeReportExpectation()))))
  
  t2 <- m1$expectation$output$thresholds

  t2 <- t2[,colnames(t1),drop=FALSE]

  expect_equivalent(is.na(t1), is.na(t2))
  
  mask <- !is.na(t1)
  expect_equal(t1[mask], t2[mask], 1e-9)
}

# from package countreg
qzinbinom <- function(p, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  p <- pmax(0, (p - pi)/(1 - pi))
  rval <- qnbinom(p, mu = mu, size = theta, lower.tail = TRUE, log.p = FALSE)
  rval
}

test_that("Normal", {
  library(OpenMx)
  library(testthat)
  
  RNGversion("4.0")
  set.seed(1)
  
  factorModel <- mxModel(
    "One Factor",
    mxMatrix(nrow=1, ncol=5, free=FALSE, values=0, name="M"),
    mxMatrix("Full", 5, 1, values=0.8,
             free=TRUE, name="A"),
    mxMatrix("Symm", 1, 1, values=1,
             free=FALSE, name="L"),
    mxMatrix("Diag", 5, 5, values=1,
             free=FALSE, name="U"),
    mxAlgebra(A %*% L %*% t(A) + U, name="R"),
    mxMatrix(nrow=3, ncol=3,
             values=c(.01, 1.1, NA,
                      .01, 2, NA,
                      .01, 4, .5),
             dimnames=list(c(), paste0('x',1:3)),
             name="D"),
    mxExpectationNormal(covariance = "R",
                        dimnames = paste0('x',1:5),
                        discrete = "D",
                        discreteSpec = matrix(c(4, 1, 5, 1, 6, 2), ncol=3,
                                              dimnames=list(c(), paste0('x',1:3))),
                        means = "M"),
    mxFitFunctionML())
  
  thr <- mxGetExpected(factorModel, "thresholds")
  expect_equal(thr[1:4,'x1'], c(-0.53, 0.68, 1.65, 2.5), .01)
  expect_equal(thr[1:5,'x2'], c(-1.36, -0.28, 0.6, 1.38, 2.08), .01)
  expect_equal(thr[1:6,'x3'], c(-1.87, -1.1, -0.49, 0.02, 0.46, 0.86), .01)
  
  factorModel <- mxGenerateData(factorModel, 400, returnModel = TRUE)

  # raw counts as integer
  factorModel$data$observed$x1 <-
    unclass(factorModel$data$observed$x1) - 1L
  
  # raw counts as numeric
  factorModel$data$observed$x2 <-
    unclass(factorModel$data$observed$x2) - 1.0

  factorModel$D$free <- !is.na(factorModel$D$values)
  
  fit <- mxRun(factorModel)
  expect_equal(fit$output$fit, 6256.76, .01)
  dv <- fit$D$values
  dv <- dv[!is.na(dv)]
  expect_equal(dv, c(0, 1.047, 0.022, 2.072, 0.013, 2.91, 0.411), .01)
  
  factorModel$expectation$discreteSpec <-
    factorModel$expectation$discreteSpec[,c(3,1,2)]
  expect_error(mxRun(factorModel),
               "must have the same column names")
  expect_error(mxGetExpected(factorModel, "thresholds"),
               "must have the same column names")
})

# ---------------

test_that("RAM", {
  library(OpenMx)
  library(testthat)
  
  RNGversion("4.0")
  set.seed(1)
  
  manifests <- paste0('x',1:5)
  latents <- c("G")
  factorModel <- mxModel(
    "One Factor",
    type="RAM",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests,values=0.8),
    mxPath(from=manifests, arrows=2,values=1, free=FALSE),
    mxPath(from=latents, arrows=2,
           free=FALSE, values=1.0),
    mxPath(from = 'one', to = manifests, free=FALSE),
    mxMarginalPoisson(paste0("x",1:2), c(4,5), c(1.1, 2)),
    mxMarginalNegativeBinomial("x3", 6, 4, .5))
  
  trueDv <- factorModel$Discrete$values
  trueDv <- trueDv[!is.na(trueDv)]
  
  thr <- mxGetExpected(factorModel, "thresholds")
  expect_equal(thr[1:4,'x1'], c(-0.53, 0.68, 1.65, 2.5), .01)
  expect_equal(thr[1:5,'x2'], c(-1.36, -0.28, 0.6, 1.38, 2.08), .01)
  expect_equal(thr[1:6,'x3'], c(-1.87, -1.1, -0.49, 0.02, 0.46, 0.86), .01)
  
  factorModel <- mxGenerateData(factorModel, 400, returnModel = TRUE)
  
  verifyFrontBackMatch(factorModel)
  
  # autodetect maximum count from data
  factorModel$expectation$discreteSpec[1,] <- NA
  
  fit <- mxRun(factorModel)
  #  summary(fit)
  expect_equal(fit$output$fit, 6256.76, .01)
  dv <- fit$Discrete$values
  dv <- dv[!is.na(dv)]
  expect_equal(dv, c(0, 1.047, 0.022, 2.072, 0.013, 2.91, 0.411), .01)

  expect_equal(colnames(fit$expectation$discreteSpec), paste0('x',1:3))  
  verifyFrontBackMatch(fit)
})

test_that("mediation", {
  library(OpenMx)
  library(testthat)
  library(MASS)

  RNGversion("4.0")
  set.seed(1)
  
  N<-5000
  A<-matrix(0,3,3)
  b_1_2<-.5
  b_2_3<-.5
  b_1_3<-.5
  
  size<- .5
  prob<- .25
  zif<-.3
  mu<- (1-prob)*size/prob
  
  A[2,1]<-b_1_2   #X1 to X2
  A[3,2]<-b_2_3   #X2 to X3
  A[3,1]<-b_1_3   #X1 to X3
  S<-diag(c(1,
            1-b_1_2^2, 
            1 - b_2_3^2*(1-b_1_2^2) - (b_2_3*b_1_2)^2 - b_1_3^2 -  2*b_1_2*b_2_3*b_1_3
  ), 3)
  R<-solve(diag(3)-A)%*%S%*%t(solve(diag(3)-A))
  
  z<-mvrnorm(N, rep(0,3), R, empirical=F)
  z[,2]<-qzinbinom(pnorm(z[,2], 0, 1), size=size, mu=mu, pi=zif)
  
  manifests<-paste0("x",1:3)
  colnames(z)<-manifests
  
  z.f<-z<-as.data.frame(z)
  z.f[,2]<-mxFactor(z[,2], levels=0:max(z[,2]))
  
  factorModel<-mxModel(
    name="countMed",type="RAM", manifestVars = manifests, latentVars =NULL,
    
    mxPath(from=c("x1","x2"), to="x3", arrows=1, free=T),
    mxPath(from="x1", to="x2", arrows=1, free=T),
    
    mxPath(from=manifests, arrows=2, values=1, free=c(T,F,T) ),
    mxPath(from="one", to=manifests, arrows=1, values=0, free=c(T,F,T)),
    
    mxMarginalNegativeBinomial(c("x2"), size=1, prob=.5),
    mxData(z.f, type="raw"))
  
  expect_error(mxGetExpected(factorModel, "thresholds"),
               "maximum observed count")

  fit <-mxRun(factorModel)
#  summary(fit)
  expect_equivalent(fit$Discrete$values, c(zif, size, prob), .025)
  expect_equivalent(fit$M$values, rep(0,3), .01)
  expect_equivalent(fit$S$values[fit$S$free], diag(S)[c(1,3)], .05)
  expect_equivalent(fit$A$values[fit$A$free], A[fit$A$free], .1)
  verifyFrontBackMatch(fit)
})

test_that("LISREL", {
  library(OpenMx)
  library(testthat)
  
  RNGversion("4.0")
  set.seed(1)
  
  manifests <- paste0('x',1:5)
  latents <- c("G")
  factorModel <- mxModel(
    "One Factor",
    type="LISREL",
    manifestVars = manifests,
    latentVars = latents,
    mxPath(from=latents, to=manifests,values=0.8),
    mxPath(from=manifests, arrows=2,values=1, free=FALSE),
    mxPath(from=latents, arrows=2,
           free=FALSE, values=1.0),
    mxPath(from = 'one', to = manifests, free=FALSE),
    mxMarginalPoisson(paste0("x",1:2), c(4,5), c(1.1, 2)),
    mxMarginalNegativeBinomial("x3", 6, 4, .5))
  
  trueDv <- factorModel$Discrete$values
  trueDv <- trueDv[!is.na(trueDv)]
  
  thr <- mxGetExpected(factorModel, "thresholds")
  expect_equal(thr[1:4,'x1'], c(-0.53, 0.68, 1.65, 2.5), .01)
  expect_equal(thr[1:5,'x2'], c(-1.36, -0.28, 0.6, 1.38, 2.08), .01)
  expect_equal(thr[1:6,'x3'], c(-1.87, -1.1, -0.49, 0.02, 0.46, 0.86), .01)
  
  factorModel <- mxGenerateData(factorModel, 400, returnModel = TRUE)
  
  verifyFrontBackMatch(factorModel)
  
  # autodetect maximum count from data
  factorModel$expectation$discreteSpec[1,] <- NA
  
  fit <- mxRun(factorModel)
  #  summary(fit)
  expect_equal(fit$output$fit, 6256.76, .01)
  dv <- fit$Discrete$values
  dv <- dv[!is.na(dv)]
  expect_equal(dv, c(0, 1.047, 0.022, 2.072, 0.013, 2.91, 0.411), .01)
  
  expect_equal(colnames(fit$expectation$discreteSpec), paste0('x',1:3))  
  verifyFrontBackMatch(fit)
})

test_that("probit+poisson ML+WLS", {
  library(OpenMx)
  library(testthat)
  
  RNGversion("4.0")
  set.seed(1)
  
  data("jointdata", package ="OpenMx")

  jointdata[,c(2,4,5)] <-
    mxFactor(jointdata[,c(2,4,5)], 
             levels=list(c(0,1), c(0, 1, 2, 3), c(0, 1, 2)))

  jointdata$extra <- rnorm(nrow(jointdata))
  
  build <- function(wls=FALSE) {
    m1 <- mxModel(
      "m1", type="RAM",
      manifestVars = paste0('z',sample.int(5,5)),
      latentVars='G',
      mxData(jointdata[,sample.int(5,5)], "raw", verbose=0L),
      mxPath('one', paste0('z', c(1,3))),
      mxPath(paste0('z', c(1,3)), arrows=2, values=2),
      mxPath(paste0('z', c(2,4,5)), arrows=2, free=FALSE, values=.5),
      mxPath('G', arrows=2, values=1, free=FALSE),
      mxPath('G', paste0('z', 1:5), values=1),
      mxMarginalProbit(paste0('z', c(2,5)), nThresh=c(1,2), free=TRUE),
      mxMarginalPoisson('z4', lambda = .5))
    if (wls) m1 <- mxModel(m1, mxFitFunctionWLS())
    m1
  }
  
  m1 <- mxRun(build())
  verifyFrontBackMatch(m1)
  m2 <- mxRun(build())
  verifyFrontBackMatch(m2)
  m3 <- mxRun(build())
  verifyFrontBackMatch(m3)
  expect_equal(m1$output$fit - m2$output$fit, 0, 1e-9)
  expect_equal(m1$output$fit - m3$output$fit, 0, 1e-9)

  m1 <- mxRun(build(TRUE))
  verifyFrontBackMatch(m1)
  m2 <- mxRun(build(TRUE))
  verifyFrontBackMatch(m2)
  m3 <- mxRun(build(TRUE))
  verifyFrontBackMatch(m3)
  expect_equal(m1$output$fit - m2$output$fit, 0, 1e-9)
  expect_equal(m1$output$fit - m3$output$fit, 0, 1e-9)
})
