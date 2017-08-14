library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 30
i1 <- rpf.drm(multidimensional=TRUE)
items <- list()
items[1:numItems] <- list(i1)
correct <- matrix(NA, 4, numItems)
for (x in 1:numItems) correct[,x] <- rpf.rparam(i1, version=1)
correct[1,] <- 1
correct[3,] <- logit(0)
correct[4,] <- logit(1)

ip.mat <- mxMatrix(name="item", nrow=4, ncol=numItems,
                   values=c(1,0, logit(0), logit(1)),
                   free=c(FALSE, TRUE, FALSE, FALSE))
colnames(ip.mat) <- paste("i", 1:numItems, sep="")
rownames(ip.mat) <- c('f1', 'b','g','u')

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=TRUE)
rownames(m.mat) <- 'f1'
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=TRUE)
dimnames(cov.mat) <- list('f1','f1')

latent <- mxModel('latent', m.mat, cov.mat,
		  mxDataDynamic("cov", expectation="latentTest.expectation"),
		  mxExpectationNormal(covariance="cov", means="mean"),
		  mxFitFunctionML())

# We generate data in the "wrong" order to preserve compatibility
# with the previous version of the test.
data <- rpf.sample(500, items, correct, cov=matrix(5,1,1))

ldata <- rpf.sample(300, items, correct, mean=.5, cov=matrix(5,1,1))

ip.fix <- ip.mat
ip.fix$free[,] <- FALSE
ip.fix$values[,] <- correct

if (0) {
  plan <- mxComputeEM('expectation', 'scores',
                      mxComputeGradientDescent(paste('latent',c('mean','cov'), sep="."),
                                               fitfunction="latent.fitfunction"),
                      verbose=0L)
}

plan <- mxComputeIterate(list(
    mxComputeOnce('expectation', 'scores'),
  mxComputeGradientDescent(paste('latent',c('mean','cov'), sep="."),
                           fitfunction="latent.fitfunction"),
    mxComputeOnce('expectation'),
  mxComputeOnce('fitfunction', 'fit')))

m1 <- mxModel(model="latentTest", ip.fix, latent,
              mxData(observed=ldata, type="raw"),
              mxExpectationBA81(items, mean="latent.mean", cov="latent.cov"),
              mxFitFunctionML(),
              plan)
m1Fit <- mxRun(m1)
omxCheckCloseEnough(m1Fit$output$estimate, c(.46, 4.44), .1)
omxCheckCloseEnough(m1Fit$output$fit, 8251.09, .01)

latent.plan <- NULL
if (1) {
	latent.plan <- mxComputeGradientDescent('latent.cov', fitfunction="latent.fitfunction")
} else {
	latent.plan <- 
	    mxComputeSequence(list(mxComputeOnce('expectation'),
				   mxComputeOnce('expectation', "latentDistribution", "copy"),
				   mxComputeOnce('fitfunction', "set-starting-values")),
			      freeSet='latent.cov')
}

latent$mean$free[,] <- FALSE
latent$data$expectation <- "drmmg.expectation"

m2 <- mxModel(model="drmmg", ip.mat, latent,
	      mxData(observed=data, type="raw"),
	      mxExpectationBA81(items, mean="latent.mean", cov="latent.cov"),
	      mxFitFunctionML(),
	      mxComputeEM('expectation', 'scores',
			  mxComputeSequence(list(
			      mxComputeNewtonRaphson(freeSet='item'),
			      latent.plan)),
			  verbose=0L))

m2 <- mxRun(m2)

omxCheckCloseEnough(m2$output$fit, 14129.94, .01)
omxCheckCloseEnough(m2$submodels$latent$matrices$cov$values[1,1], 4.377, .01)
		
emstat <- m2$compute$output
omxCheckCloseEnough(emstat$EMcycles, 29, 4)
#omxCheckCloseEnough(emstat$totalMstep, 763, 20) # includes the latent distribution

					#print(m2$matrices$item$values)
					#print(correct.mat)
mask <- is.finite(correct)
got <- cor(c(m2$matrices$item$values[mask]),
	   c(correct[mask]))
omxCheckCloseEnough(got, .994, .01)

if (1) {
  ip.mat <- mxMatrix(name="item", nrow=4, ncol=numItems,
                     values=c(1,0, logit(0), logit(1)),
                     free=c(TRUE, TRUE, FALSE, FALSE))
  colnames(ip.mat) <- paste("i", 1:numItems, sep="")
  rownames(ip.mat) <- c('f1', 'b','g','u')
  ip.mat$labels[1,] <- 'a1'
  
  m2 <- mxModel(model="drmmg", ip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=items),
                mxFitFunctionML(),
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation', 'scores'),
                  mxComputeOnce('fitfunction', c('gradient', 'hessian', 'ihessian')),
                  mxComputeReportDeriv()
                )))
  deriv <- mxRun(m2, silent=TRUE)
  omxCheckCloseEnough(deriv$output$ihessian, solve(deriv$output$hessian), 1e-4)
  
  if (0) {
    m3 <- mxModel(m2, mxComputeSequence(list(
      mxComputeOnce('fitfunction', 'fit'),
      mxComputeNumericDeriv(parallel=FALSE, iterations=2L),
      mxComputeReportDeriv())))
    deriv <- mxRun(m3)
    stop("ok")
  }
  
  m2 <- mxModel(model="drmmg", ip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(ItemSpec=items),
                mxFitFunctionML(),
                mxComputeEM('expectation', 'scores',
                            mxComputeNewtonRaphson(freeSet='item')))

  m2 <- mxRun(m2)
  emstat <- m2$compute$output
  omxCheckCloseEnough(emstat$EMcycles, 21, 1)
  omxCheckCloseEnough(emstat$totalMstep, 56, 5)
  omxCheckCloseEnough(m2$fitfunction$result, 14129.04, .01)
  omxCheckCloseEnough(m2$matrices$item$values[1,], rep(2.133, numItems), .002)
  # correct values are from flexMIRT
  est <- c(-0.838622, -1.02653, -0.0868472, -0.251784, 0.953364,  0.735258, 0.606918,
           1.04239, 0.466055, -2.05196, -0.0456446,  -0.320668, -0.362073, 2.02502,
           0.635298, -0.0731132, -2.05196,  -0.0456446, -1.17429, 0.880002, -0.838622,
           -0.838622, 1.02747,  0.424094, -0.584298, 0.663755, 0.663755, 0.064287, 1.38009,
           1.01259 )
  omxCheckCloseEnough(m2$matrices$item$values[2,], est, .002)
}

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
  
  # MH-RM takes forever, not run
  pars <- confmirt(rdata, 1, itemtype="2PL", D=1, quadpts=49, pars='values')
  pars[pars$name=="a1",'value'] <- 1
  pars[pars$name=="a1",'est'] <- FALSE
  pars[pars$name=="COV_11",'est'] <- TRUE
  fit <- confmirt(rdata, 1, itemtype="2PL", D=1, quadpts=49, pars=pars)
  got <- coef(fit)
  got$GroupPars <- NULL
  round(m2$matrices$item$values - sapply(got, function(l) l[1,]), 2)
}
