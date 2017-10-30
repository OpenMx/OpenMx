#options(error = browser)
library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 20
items <- list()
items[1:numItems] <- list(rpf.grm(outcomes=2))  # equivalent to 2PL

# create random item parameters
correct <- sapply(items, rpf.rparam, version=1)

# simulate data
data.g1 <- rpf.sample(500, items, correct)
data.g2 <- rpf.sample(500, items, correct, mean=-1, cov=matrix(5,1,1))
data.g3 <- rpf.sample(500, items, correct, mean=1, cov=matrix(.5,1,1))

if (0) {
  # for flexMIRT, write CSV
  write.table(sapply(data.g1, unclass)-1, file="drm-mg2-g1.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g2, unclass)-1, file="drm-mg2-g2.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(sapply(data.g3, unclass)-1, file="drm-mg2-g3.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# This create an IFA model for each group. Item parameters are
# constrained equal across groups.
mkgroup <- function(model.name, data, latent.free) {
  ip.mat <- mxMatrix(name="item", nrow=2, ncol=numItems,
                     values=c(1,0), free=TRUE)
  colnames(ip.mat) <- colnames(data)
  rownames(ip.mat) <- c('f1', 'b')
  
  for (ix in 1:numItems) {
    for (px in 1:2) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat$labels[px,ix] <- name
    }
  }
  
  dims <- 1
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=dims, values=0)
  rownames(m.mat) <- "f1"
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims))
  dimnames(cov.mat) <- list("f1", "f1")

  mean <- "mean"
  cov <- "cov"
  if (latent.free) {
    lm <- paste(model.name, "latent", sep="")
    mean <- paste(lm, "expMean", sep=".")
    cov <- paste(lm, "expCov", sep=".")
  }
  
  m1 <- mxModel(model=model.name, ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=items, mean=mean, cov=cov,
                                  verbose=ifelse(latent.free, 0L, 0L)),
                mxFitFunctionML())
  m1
}

g1 <- mkgroup("g1", data.g1, FALSE)
g2 <- mkgroup("g2", data.g2, TRUE)
g3 <- mkgroup("g3", data.g3, TRUE)

groups <- paste("g", 1:3, sep="")

# Cannot test derivatives at starting values because Hessian starts very close to singular.

# This create a latent distribution model that can be used to impose
# equality constraints on latent distribution parameters.
mklatent <- function(name) {
	mMat <- mxMatrix(nrow=1, ncol=1, free=T, values=0, name="expMean")
	rownames(mMat) <- "f1"
	covMat <- mxMatrix(type="Symm", nrow=1, ncol=1, free=T, values=1, name="expCov")
	dimnames(covMat) <- list("f1", "f1")

	m1 <- mxModel(paste(name, "latent", sep=""), mMat, covMat,
		      mxDataDynamic("cov", expectation=paste(name, "expectation", sep=".")),
		      mxExpectationNormal(covariance="expCov", means="expMean"),
		      mxFitFunctionML())
  m1
}

latent <- mxModel("latent",
                  mxFitFunctionMultigroup(paste(paste(groups[-1],"latent",sep=""), "fitfunction", sep=".")))

g2.latent <- mklatent("g2")
g3.latent <- mklatent("g3")

latent.vargroup <- apply(expand.grid(paste(groups[-1], "latent", sep=""), c('expMean','expCov')),
                         1, paste, collapse='.')

latent.plan <- NULL  # need a plan for latent distribution parameters

if (0) {
  # Copy latent distribution parameters from current estimates without transformation.
  latent.plan <- mxComputeSequence(list(mxComputeOnce(paste(groups, 'expectation', sep='.')),
						      mxComputeOnce(paste(groups, 'expectation', sep='.'),
                                                      "latentDistribution", "copy"),  # c('mean','covariance')
                                        mxComputeOnce('fitfunction', "set-starting-values")),
                                   freeSet=latent.vargroup)
} else {
  # Obtain latent distribution parameters via mxExpectationNormal.
  # This permits equality constraints (and potentially more complex latent structure).
#  latent.plan <- mxComputeGradientDescent(latent.vargroup, fitfunction="latent.fitfunction")
  latent.plan <- mxComputeNewtonRaphson(latent.vargroup, fitfunction="latent.fitfunction")
}

grpModel <- mxModel(model="groupModel", g1, g2, g3, g2.latent, g3.latent, latent,
                    mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                    mxComputeSequence(list(
                      mxComputeEM(paste(groups, 'expectation', sep='.'), 'scores',
                                  mxComputeSequence(list(
				      mxComputeNewtonRaphson(freeSet=paste(groups,'item',sep=".")),
				      latent.plan)),
                                  information="mr1991", verbose=0L,
				  tolerance=1e-10,
                                  infoArgs=list(fitfunction=c("fitfunction", "latent.fitfunction"))),
                      mxComputeStandardError(),
                      mxComputeHessianQuality(),
                    mxComputeOnce('fitfunction', 'gradient'),
                   mxComputeReportDeriv())))

  #grpModel <- mxOption(grpModel, "Number of Threads", 1)
  
grpModel <- mxOption(grpModel, "Checkpoint Units", "iterations")
grpModel <- mxOption(grpModel, "Checkpoint Count", 1)
grpModel <- mxRun(grpModel)  #, checkpoint = TRUE

#dm <- grpModel$compute$steps[[1]]$debug$rateMatrix

omxCheckCloseEnough(max(abs(grpModel$output$gradient)), 0, .17)

omxCheckCloseEnough(grpModel$output$fit, 30114.94, .02)
omxCheckEquals(grpModel$output$fitUnits, "-2lnL")
omxCheckCloseEnough(AIC(logLik(grpModel)), 30202.95, .02)
omxCheckCloseEnough(BIC(logLik(grpModel)), 30436.73, .02)
omxCheckCloseEnough(summary(grpModel)$informationCriteria['AIC:','par'], 30202.95, .02)
omxCheckCloseEnough(summary(grpModel)$informationCriteria['BIC:','par'], 30436.73, .02)

  omxCheckCloseEnough(grpModel$submodels$g2latent$matrices$expMean$values, -.834, .01)
  omxCheckCloseEnough(grpModel$submodels$g2latent$matrices$expCov$values, 3.93, .01)
  omxCheckCloseEnough(grpModel$submodels$g3latent$matrices$expMean$values, .933, .01)
  omxCheckCloseEnough(grpModel$submodels$g3latent$matrices$expCov$values, .444, .01)

emstat <- grpModel$compute$steps[[1]]$output
if (0) {
# Optimization path has too much variance
omxCheckCloseEnough(emstat$EMcycles, 50, 4)
omxCheckCloseEnough(emstat$totalMstep, 380, 12)
omxCheckCloseEnough(emstat$semProbeCount, 90, 15)
  
#  cat(deparse(round(grpModel$output$standardErrors, 3)))
  semse <- c(0.069, 0.077, 0.074, 0.077, 0.094, 0.097, 0.125,  0.111, 0.069, 0.074,
             0.132, 0.116, 0.08, 0.081, 0.209, 0.163,  0.102, 0.133, 0.114, 0.107,
             0.205, 0.151, 0.068, 0.077, 0.073,  0.138, 0.078, 0.081, 0.088, 0.087,
             0.061, 0.068, 0.125, 0.11,  0.084, 0.09, 0.094, 0.094, 0.092, 0.089,
             0.11, 0.399, 0.068,  0.055)
  omxCheckCloseEnough(c(grpModel$output$standardErrors), semse, .03)
  omxCheckCloseEnough(log(grpModel$output$conditionNumber), 5.1, .5)
omxCheckTrue(grpModel$output$infoDefinite)
}

refModels <- mxRefModels(grpModel, run=TRUE)
ind <- refModels[['Independence']]
omxCheckCloseEnough(ind$output$fit, 38588.19, .01)
omxCheckCloseEnough(summary(ind)$informationCriteria['AIC:','par'], 38628.187, .02)
omxCheckCloseEnough(summary(ind)$informationCriteria['BIC:','par'], 38734.451, .02)

i1 <- mxModel(grpModel,
                mxComputeSequence(steps=list(
                  mxComputeOnce('fitfunction', 'information', "meat"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
i1 <- mxRun(i1)
  
  #cat(deparse(round(i1$output$standardErrors,3)))
  se <- c(0.071, 0.078, 0.076, 0.079, 0.097, 0.099, 0.132,  0.117, 0.075,
          0.077, 0.135, 0.121, 0.081, 0.083, 0.215, 0.169,  0.111, 0.141,
          0.121, 0.113, 0.213, 0.159, 0.074, 0.082, 0.077,  0.139, 0.084,
          0.087, 0.095, 0.09, 0.064, 0.07, 0.135, 0.115,  0.091, 0.095, 0.097,
          0.098, 0.096, 0.093, 0.12, 0.512, 0.072,  0.057)
  omxCheckCloseEnough(c(i1$output$standardErrors), se, .01)
  omxCheckCloseEnough(log(i1$output$conditionNumber), 5.6, .5)

if (0) {
  library(mirt)
  rdata <- sapply(data, unclass)-1
  # for flexMIRT, write CSV
  #write.table(rdata, file="ifa-drm-mg.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
  pars <- mirt(rdata, 1, itemtype="2PL", D=1, quadpts=49, pars='values')
  pars[pars$name=="a1",'value'] <- 1
  pars[pars$name=="a1",'est'] <- FALSE
  pars[pars$name=="COV_11",'est'] <- TRUE
  fit <- mirt(rdata, 1, itemtype="2PL", D=1, quadpts=49, pars=pars)
  # LL -7064.519 * -2 = 14129.04
  got <- coef(fit)
  print(got$GroupPars)
  # COV 4.551
  got$GroupPars <- NULL
  round(m2$matrices$item$values - simplify2array(got), 2)
}
