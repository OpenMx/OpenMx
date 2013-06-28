options(error = utils::recover)
library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 10
i1 <- rpf.drm(multidimensional=TRUE)
items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
  correct[[ix]][3] <- 0
  correct[[ix]][4] <- 1
}
correct.mat <- simplify2array(correct)
correct.mat[2,] <- correct.mat[2,] * -correct.mat[1,]

ability <- rnorm(500)
data <- rpf.sample(ability, items, correct.mat)

spec <- mxMatrix(name="ItemSpec", nrow=3, ncol=numItems,
         values=sapply(items, function(m) slot(m,'spec')),
         free=FALSE, byrow=TRUE)

ip.mat <- mxMatrix(name="itemParam", nrow=4, ncol=numItems,
                   values=c(1,0,0, 1),
                   free=c(TRUE, TRUE, FALSE, FALSE),
		   lbound=c(1e-6, -1e6, 0, 0))
ip.mat@free.group <- 'param'

eip.mat <- mxMatrix(name="EItemParam", nrow=4, ncol=numItems,
                    values=c(1,0,0, 1),
                    free=c(TRUE, TRUE, FALSE, FALSE))

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

m2 <- mxModel(model="drm1", ip.mat, spec, m.mat, cov.mat, eip.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec="ItemSpec", EItemParam="EItemParam",
                mean="mean", cov="cov",
		qpoints=30,
		scores="full"),
              mxFitFunctionBA81(ItemParam="itemParam"),
	      mxComputeIterate(steps=list(
				 mxComputeAssign(from="itemParam", to="EItemParam"),
				 mxComputeOnce(expectation='expectation', context='E'),
				 mxComputeGradientDescent(free.group='param'),
				 mxComputeOnce(expectation='expectation', context='M'),
				 mxComputeOnce(fitfunction='fitfunction')
				 )))

	m2 <- mxOption(m2, "Analytic Gradients", 'Yes')
	m2 <- mxOption(m2, "Verify level", '-1')
m2 <- mxOption(m2, "Function precision", '1.0E-5')
m2 <- mxRun(m2)

#print(m2@matrices$itemParam@values)
#print(correct.mat)
got <- cor(c(m2@matrices$itemParam@values[1:2,]),
           c(correct.mat[1:2,]))
omxCheckCloseEnough(got, .988, .01)
ability <- scale(ability)
omxCheckCloseEnough(m2@output$ability[1,], as.vector(ability), 3.5*max(m2@output$ability[,2]))
omxCheckCloseEnough(cor(c(m2@output$ability[1,]), ability), .737, .01)
