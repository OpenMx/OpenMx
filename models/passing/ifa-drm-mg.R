library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 30
i1 <- rpf.drm(multidimensional=TRUE)
items <- list()
items[1:numItems] <- i1
correct <- matrix(NA, 4, numItems)
for (x in 1:numItems) correct[,x] <- rpf.rparam(i1)
correct[1,] <- 1
correct[3,] <- 0
correct[4,] <- 1

data <- rpf.sample(500, items, correct, cov=matrix(5,1,1))

	ip.mat <- mxMatrix(name="itemParam", nrow=4, ncol=numItems,
			   values=c(1,0,0, 1),
			   free=c(FALSE, TRUE, FALSE, FALSE))
	
	m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
	cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=TRUE)

	m2 <- mxModel(model="drmmg", ip.mat, m.mat, cov.mat,
	              mxData(observed=data, type="raw"),
	              mxExpectationBA81(mean="mean", cov="cov",
	                                ItemSpec=items, ItemParam="itemParam"),
	              mxFitFunctionML(),
	              mxComputeEM('expectation', 'scores',
	                          mxComputeNewtonRaphson(free.set='itemParam'),
	                          mxComputeSequence(list(mxComputeOnce('expectation', "latentDistribution", "copy"),
	                                                 mxComputeOnce('fitfunction', "starting")),
	                                            free.set='cov'),
	                          mxComputeOnce('fitfunction', 'fit')))
	
	if (0) {
		fm <- read.flexmirt("/home/joshua/irt/ifa-drm-mg/ifa-drm-mg-prm.txt")
		cModel <- m2
		cModel@matrices$itemParam@values[2,] <- fm$G1$param[2,]
		cModel@matrices$cov@values <- fm$G1$cov
		cModel <- mxModel(cModel,
				  mxExpectationBA81(mean="mean", cov="cov",
						    ItemSpec="ItemSpec",
						    scores="full"),
				  mxComputeSequence(steps=list(
						      mxComputeOnce('expectation'),
						      mxComputeOnce('fitfunction', 'fit'))))
		cModel <- mxRun(cModel)
		cModel@matrices$cov@values - fm$G1$cov
		cModel@output$minimum
	}

# 		m2 <- mxOption(m2, "Analytic Gradients", 'Yes')
# 		m2 <- mxOption(m2, "Verify level", '-1')
# 		m2 <- mxOption(m2, "Function precision", '1.0E-5')
		m2 <- mxRun(m2)

emstat <- m2@compute@output
omxCheckCloseEnough(emstat$EMcycles, 12, 1)
omxCheckCloseEnough(emstat$totalMstep, 35, 2)

omxCheckCloseEnough(m2@output$fit, 14129.94, .01)
		omxCheckCloseEnough(m2@matrices$cov@values[1,1], 4.377, .01)
		
					#print(m2@matrices$itemParam@values)
					#print(correct.mat)
		got <- cor(c(m2@matrices$itemParam@values),
			   c(correct))
		omxCheckCloseEnough(got, .994, .01)

if (1) {
  ip.mat <- mxMatrix(name="itemParam", nrow=4, ncol=numItems,
                     values=c(1,0,0, 1),
                     free=c(TRUE, TRUE, FALSE, FALSE))
  ip.mat@labels[1,] <- 'a1'
  
  m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
  cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

  m2 <- mxModel(model="drmmg", ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=items, ItemParam="itemParam"),
                mxFitFunctionML(),
                mxComputeSequence(steps=list(
                  mxComputeOnce('expectation', 'scores'),
                  mxComputeOnce('fitfunction', c('gradient', 'hessian', 'ihessian'))
                )))
  m2 <- mxRun(m2)
  omxCheckCloseEnough(m2@output$ihessian, solve(m2@output$hessian), 1e-4)
  
  m2 <- mxModel(model="drmmg", ip.mat, m.mat, cov.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
                                  ItemSpec=items, ItemParam="itemParam"),
                mxFitFunctionML(),
                mxComputeEM('expectation', 'scores',
                            mxComputeNewtonRaphson(free.set='itemParam'),
                            mxComputeNothing(),
                            mxComputeOnce('fitfunction', 'fit')))

  m2 <- mxRun(m2)
  emstat <- m2@compute@output
  omxCheckCloseEnough(emstat$EMcycles, 38, 1)
  omxCheckCloseEnough(emstat$totalMstep, 503, 10)
  omxCheckCloseEnough(m2@fitfunction@result, 14129.04, .01)
  omxCheckCloseEnough(m2@matrices$itemParam@values[1,], rep(2.133, numItems), .002)
  # correct values are from flexMIRT
  est <- c(-0.838622, -1.02653, -0.0868472, -0.251784, 0.953364,  0.735258, 0.606918,
           1.04239, 0.466055, -2.05196, -0.0456446,  -0.320668, -0.362073, 2.02502,
           0.635298, -0.0731132, -2.05196,  -0.0456446, -1.17429, 0.880002, -0.838622,
           -0.838622, 1.02747,  0.424094, -0.584298, 0.663755, 0.663755, 0.064287, 1.38009,
           1.01259 )
  omxCheckCloseEnough(m2@matrices$itemParam@values[2,], est, .002)
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
  round(m2@matrices$itemParam@values - simplify2array(got), 2)
  
  # MH-RM takes forever, not run
  pars <- confmirt(rdata, 1, itemtype="2PL", D=1, quadpts=49, pars='values')
  pars[pars$name=="a1",'value'] <- 1
  pars[pars$name=="a1",'est'] <- FALSE
  pars[pars$name=="COV_11",'est'] <- TRUE
  fit <- confmirt(rdata, 1, itemtype="2PL", D=1, quadpts=49, pars=pars)
  got <- coef(fit)
  got$GroupPars <- NULL
  round(m2@matrices$itemParam@values - sapply(got, function(l) l[1,]), 2)
}
