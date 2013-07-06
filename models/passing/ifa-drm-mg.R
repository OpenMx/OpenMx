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

if (1) {
  spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
                   values=sapply(items, function(m) slot(m,'spec')),
                   free=FALSE, byrow=TRUE)
  
  ip.mat <- mxMatrix(name="itemParam", nrow=4, ncol=numItems,
                     values=c(1,0,0, 1),
                     free=c(FALSE, TRUE, FALSE, FALSE))
  ip.mat@free.group <- 'param'
  
  eip.mat <- mxMatrix(name="EItemParam", nrow=4, ncol=numItems,
		      values=c(1,0,0, 1),
		      free=TRUE)

  m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
  cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=TRUE)

  m2 <- mxModel(model="drmmg", ip.mat, spec, m.mat, cov.mat, eip.mat,
                mxData(observed=data, type="raw"),
                mxExpectationBA81(mean="mean", cov="cov",
				  ItemSpec="ItemSpec",
				  EItemParam="EItemParam"),
                mxFitFunctionBA81(ItemParam="itemParam"),
		# integrate FitFuncBA81 into FitFuncML
		mxComputeIterate(steps=list(
				   mxComputeAssign(from="itemParam", to="EItemParam"), # replace with "fixed" mxAlgebra
				   mxComputeOnce('expectation', context='E'),
				   mxComputeGradientDescent(free.group='param'),
				   # list=(matrix, free parameter, models)
				   # can recompute just from dependencies?
				   mxComputeOnce('expectation', context='M'),
				   mxComputeOnce('fitfunction')
				 )))
  
  m2 <- mxOption(m2, "Analytic Gradients", 'Yes')
  m2 <- mxOption(m2, "Verify level", '-1')
  m2 <- mxOption(m2, "Function precision", '1.0E-7')
  m2 <- mxRun(m2)
  
#  omxCheckCloseEnough(m2@output$minimum, 14130.2, 1) TODO
  omxCheckCloseEnough(m2@matrices$cov@values[1,1], 4.377, .01)
  
  #print(m2@matrices$itemParam@values)
  #print(correct.mat)
  got <- cor(c(m2@matrices$itemParam@values),
             c(correct))
  omxCheckCloseEnough(got, .994, .01)
}

if (0) {
  library(mirt)
  rdata <- sapply(data, unclass)-1
  # for flexMIRT, write CSV
  write.table(rdata, file="irt-drm-mg.csv", quote=FALSE, row.names=FALSE, col.names=FALSE)
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
