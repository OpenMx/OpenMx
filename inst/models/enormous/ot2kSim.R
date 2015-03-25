library(ggplot2)
library(OpenMx)
library(rpf)
library(pROC)

mkmodel <- function(numItems, numPeople, genCond, fitCond, numBad) {
  spec <- list()
  
  trueFit <- rep(FALSE, numBad)
  trueFit <- c(trueFit, rep(TRUE, numItems - length(trueFit)))
  
  if (genCond == "2pl") {
    spec[1:numItems] <- rpf.grm()
    correct <- sapply(spec, rpf.rparam)
    correct['a', trueFit==TRUE] <- 1
    data <- rpf.sample(numPeople, spec, correct)
    
    if (fitCond == "1pl") {
      spec[1:numItems] <- rpf.grm()
      ip.mat <- mxMatrix(name="ItemParam", nrow=2, ncol=numItems,
                         values=c(1,0),
                         free=c(FALSE, TRUE),
                         dimnames=list(rownames(correct), colnames(data)))
    } else { stop(fitCond) }
  }
  
  if (genCond == "3pl") {
    spec[1:numItems] <- rpf.drm()
    correct <- sapply(spec, rpf.rparam)
    correct['g',trueFit==TRUE] <- logit(0)
    correct['u',trueFit==TRUE] <- logit(1)
    correct['u',trueFit==FALSE] <- logit(.5)
    correct['u',trueFit==FALSE] <- logit(1)
    data <- rpf.sample(numPeople, spec, correct)

    if (fitCond == "2pl") {
      spec[1:numItems] <- rpf.drm()
      ip.mat <- mxMatrix(name="ItemParam", nrow=4, ncol=numItems,
                         values=c(1,0, logit(0), logit(1)),
                         free=c(TRUE, TRUE, FALSE, FALSE),
                         dimnames=list(rownames(correct), colnames(data)))
    } else { stop(fitCond) }
  }
  
  if (genCond == "nom2") {
    spec[1:numItems] <- rpf.nrm(outcomes=4)
    correct <- sapply(spec, rpf.rparam)
    half <- numItems/2
    correct[c('alf2', 'alf3'), trueFit==TRUE] <- 0
    data <- rpf.sample(numPeople, spec, correct)
    
    if (fitCond == "grm") {
      spec[1:numItems] <- rpf.grm(outcomes=4)
      ip.mat <- mxMatrix(name="ItemParam", nrow=4, ncol=numItems,
                         values=rpf.rparam(spec[[1]]),
                         free=TRUE,
                         dimnames=list(names(rpf.rparam(spec[[1]])), colnames(data)))
    } else { stop(fitCond) }
  }
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE, dimnames=list("a","a"))
  cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE, dimnames=list("a","a"))
  
  m1 <- mxModel(model="ot2k", ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=spec, ItemParam="ItemParam",
                                  mean="mean", cov="cov"),
                mxFitFunctionML(),
                mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson()))
  list(model=m1, trueFit=trueFit)
}

trial <- function(numItems, numPeople, genCond, fitCond, numBad) {
  got <- mkmodel(numItems, numPeople, genCond, fitCond, numBad)
  model <- got$model
  trueFit <- got$trueFit
  model <- mxRun(model, silent=TRUE)
  
  grp <- list(spec=model$expectation$ItemSpec,
              param=model$ItemParam$values,
              mean=model$mean$values,
              cov=model$cov$values,
              data=model$data$observed)
  
  result <- expand.grid(method=c("pearson", "rms"), alt=c(TRUE, FALSE), item=1:numItems, trueFit=NA, pval=NA)
  result <- result[!(result$method=="rms" & result$alt),]
  for (ix in 1:numItems) {
    result[result$item == ix, 'trueFit'] <- trueFit[ix]
  }

  got <- rpf.SitemFit(grp, method="pearson")
  mask <- result$method=="pearson" & !result$alt
  result[mask,'pval'] <- sapply(got, function(x) x$pval)
  
  got <- rpf.SitemFit(grp, method="pearson", alt=TRUE)
  mask <- result$method=="pearson" & result$alt
  result[mask,'pval'] <- sapply(got, function(x) x$pval)
  
  got <- rpf.SitemFit(grp, method="rms")
  mask <- result$method=="rms" & !result$alt
  result[mask,'pval'] <- sapply(got, function(x) x$pval)
  result
}

result <- NULL
for (replication in 1:20) {
  result <- rbind(result, trial(20, 500, "nom2", "grm", 20))
  result <- rbind(result, trial(20, 500, "nom2", "grm", 0))
#  result <- rbind(result, trial(10, 1000, "3pl", "2pl", 1))
}

alphaPlot <- function(result) {
  aResult <- expand.grid(method=unique(result$method), alt=unique(result$alt), alpha=seq(.01,.1, .005), got=0)
  aResult$label <- paste(aResult$method, aResult$alt, sep="+")
  # need at least 10/min(alpha) replications for good accuracy
  minRows <- 10/min(aResult$alpha)
  for (r in 1:nrow(aResult)) {
    mre <- subset(result, method==aResult$method[r] & alt==aResult$alt[r] & trueFit)
    if (aResult$alpha[r] == min(aResult$alpha) &&
          nrow(mre) && nrow(mre) < minRows) {
      warning(paste("Only", nrow(mre),"rows for", aResult$label[r], "need", minRows))
    }
    aResult$got[r] <- sum(mre$pval < aResult$alpha[r]) / nrow(mre)
  }
  aResult <- aResult[is.finite(aResult$got),]
  ggplot(aResult, aes(alpha, got, color=label)) + geom_line() +
    geom_abline(intercept=0, slope=1, color="green") +
    labs(x="expected alpha", y="empirical alpha") + ylim(0,max(aResult$got))
}

if (0) {
  numItems = 20
  numPeople = 1000
  genCond = "3pl"
  fitCond = "2pl"
  roc(trueFit ~ pval, subset(result, method=="pearson" & !alt), plot=TRUE, ci=TRUE)
  roc(trueFit ~ pval, subset(result, method=="rms" & !alt), plot=TRUE, ci=TRUE)
  roc(trueFit ~ pval, subset(result, method=="pearson" & alt), plot=TRUE, ci=TRUE)  # difference at small # of items
  alphaPlot(result)
}

# bad items make other items bad so there are two feasible ways to compose the items,
# 100% bad, 100% good
# 1 bad and the rest good

# rms has no advantage in all good vs all bad, dichotomous

# rms has less power but better alpha level (who cares?)

# pearson & alt=TRUE has slightly less power than pearson & alt=FALSE with less than 10 items
