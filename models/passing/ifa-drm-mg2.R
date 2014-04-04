#options(error = browser)
library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 20
items <- list()
items[1:numItems] <- rpf.grm(outcomes=2)  # equivalent to 2PL

# create random item parameters
correct <- sapply(items, rpf.rparam)

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
  ip.mat <- mxMatrix(name="ItemParam", nrow=2, ncol=numItems,
                     values=c(1,0), free=TRUE)
  
  for (ix in 1:numItems) {
    for (px in 1:2) {
      name <- paste(c('p',ix,',',px), collapse='')
      ip.mat$labels[px,ix] <- name
    }
  }
  
  dims <- 1
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=dims, values=0)
  cov.mat <- mxMatrix(name="cov", nrow=dims, ncol=dims, values=diag(dims))

  mean <- "mean"
  cov <- "cov"
  if (latent.free) {
    lm <- paste(model.name, "latent", sep="")
    mean <- paste(lm, "expMean", sep=".")
    cov <- paste(lm, "expCov", sep=".")
  }
  
  m1 <- mxModel(model=model.name, ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=items[[1]], ItemParam="ItemParam", mean=mean, cov=cov,
                                  verbose=ifelse(latent.free, 0L, 0L)),
                mxFitFunctionML())
  m1
}

g1 <- mkgroup("g1", data.g1, FALSE)
g2 <- mkgroup("g2", data.g2, TRUE)
g3 <- mkgroup("g3", data.g3, TRUE)

groups <- paste("g", 1:3, sep="")

# Cannot test derivatives at starting values because Hessian starts very close to singular.

if(0) {
  # for S-EM debugging
  plan <- mxComputeEM(paste(groups, 'expectation', sep='.'),
                      mxComputeNewtonRaphson(free.set=paste(groups,'ItemParam',sep=".")),
                      mxComputeOnce('fitfunction', 'fit',
                                    free.set=apply(expand.grid(groups, c('mean','cov')), 1, paste, collapse='.')),
                      information=TRUE, info.method="meat", semDebug=TRUE, semMethod=seq(.001, .02, length.out=30))
}

# This create a latent distribution model that can be used to impose
# equality constraints on latent distribution parameters.
mklatent <- function(name) {
	m1 <- mxModel(paste(name, "latent", sep=""),
		      mxMatrix(nrow=1, ncol=1, free=T, values=0, name="expMean"),
		      mxMatrix(type="Symm", nrow=1, ncol=1, free=T, values=1, name="expCov"),
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
  latent.plan <- mxComputeSequence(list(mxComputeOnce(paste(groups, 'expectation', sep='.'),
                                                      "latentDistribution", "copy"),  # c('mean','covariance')
                                        mxComputeOnce('fitfunction', "set-starting-values")),
                                   free.set=latent.vargroup)
} else {
  # Obtain latent distribution parameters via mxExpectationNormal.
  # This permits equality constraints (and potentially more complex latent structure).
  latent.plan <- mxComputeGradientDescent(latent.vargroup, fitfunction="latent.fitfunction")
}

grpModel <- mxModel(model="groupModel", g1, g2, g3, g2.latent, g3.latent, latent,
                    mxFitFunctionMultigroup(paste(groups, "fitfunction", sep=".")),
                    mxComputeSequence(list(
                      mxComputeEM(paste(groups, 'expectation', sep='.'), 'scores',
                                  mxComputeNewtonRaphson(free.set=paste(groups,'ItemParam',sep=".")),
                                  latent.plan,
                                  mxComputeOnce('fitfunction', 'fit'),
                                  information=TRUE, tolerance=1e-5, verbose=0L,
				  infoArgs=list(fitfunction=c("fitfunction", "latent.fitfunction"))),
                      mxComputeStandardError(),
                      mxComputeHessianQuality())))

  #grpModel <- mxOption(grpModel, "Number of Threads", 1)
  
grpModel <- mxRun(grpModel)

#dm <- grpModel$compute$steps[[1]]$debug$rateMatrix

plot_em_map <- function(model, cem) {   # for S-EM debugging
  require(ggplot2)
  phl <- cem$debug$paramHistLen
  probeOffset <- cem$debug$probeOffset
  semDiff <- cem$debug$semDiff

  modelfit <- NULL
  result <- data.frame()
  for (vx in 1:length(model$output$estimate)) {
    len <- phl[vx]
    offset <- probeOffset[1:len, vx]
    dd <- semDiff[1:(len-1), vx]
    mid <- offset[1:(len-1)] + diff(offset)/2
    upper <- 20
    mask <- abs(diff(offset)) < .01 & dd < upper
    df <- data.frame(mid=mid[mask], diff=dd[mask])
    m1 <- lm(diff ~ 1 + I(1/mid^2), data=df)
    modelfit <- c(modelfit, summary(m1)$r.squ)
    df$model <- predict(m1)
    result <- rbind(result, cbind(vx=vx, vname=names(model$output$estimate)[vx], df))
  }
  print(mean(modelfit))
  ggplot(subset(result, vx %in% order(modelfit)[1:9])) +
    geom_point(aes(mid, diff), size=2) + geom_line(aes(mid, model), color="green") +
    facet_wrap(~vname) + labs(x="x midpoint") + ylim(0,5)
}

if (0) {
  plot_em_map(grpModel, grpModel$compute)
}

omxCheckCloseEnough(grpModel$output$fit, 30114.94, .02)
  omxCheckCloseEnough(grpModel$submodels$g2latent$matrices$expMean$values, -.834, .01)
  omxCheckCloseEnough(grpModel$submodels$g2latent$matrices$expCov$values, 3.93, .01)
  omxCheckCloseEnough(grpModel$submodels$g3latent$matrices$expMean$values, .933, .01)
  omxCheckCloseEnough(grpModel$submodels$g3latent$matrices$expCov$values, .444, .01)

emstat <- grpModel$compute$steps[[1]]$output
if (0) {
# Optimization path has too much variance
omxCheckCloseEnough(emstat$EMcycles, 127, 3)
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

i1 <- mxModel(grpModel,
                mxComputeSequence(steps=list(
                  mxComputeOnce(paste(groups, 'expectation', sep='.')),
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
  omxCheckCloseEnough(log(i1$output$conditionNumber), 5.6, .2)

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
  round(m2$matrices$itemParam$values - simplify2array(got), 2)
}
