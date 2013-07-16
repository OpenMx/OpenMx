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

ip.mat <- mxMatrix(name="itemParam", nrow=4, ncol=numItems,
                   values=c(1,0,0, 1),
                   free=c(TRUE, TRUE, FALSE, FALSE))

eip.mat <- mxAlgebra(itemParam, name="EItemParam")

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)

m2 <- mxModel(model="drm1", ip.mat, m.mat, cov.mat, eip.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec=items, EItemParam="EItemParam",
                mean="mean", cov="cov",
		qpoints=31,
		scores="full"),
              mxFitFunctionBA81(ItemParam="itemParam"),
	      mxComputeIterate(steps=list(
				 mxComputeOnce("EItemParam"),
				 mxComputeOnce('expectation', context='EM'),
				 mxComputeNewtonRaphson(free.set='itemParam'),
				 mxComputeOnce('expectation'),
				 mxComputeOnce('fitfunction')
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
scores <- m2@expectation@scores.out
omxCheckCloseEnough(scores[,1], as.vector(ability), 3.5*max(scores[,2]))
omxCheckCloseEnough(cor(c(scores[,1]), ability), .737, .01)
