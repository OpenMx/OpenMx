#options(error = utils::recover)
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
                ItemSpec=items, EItemParam="EItemParam",ItemParam="itemParam",
                mean="mean", cov="cov", qpoints=31),
              mxFitFunctionML(),
	      mxComputeOnce('expectation', context='EM'))
m2 <- mxRun(m2)
omxCheckCloseEnough(fivenum(m2@expectation@patternLikelihood),
                    c(0.0005285109, 0.0006143148, 0.0006143148, 0.0009690677, 0.0289785643), .001)
omxCheckCloseEnough(sum(m2@expectation@patternLikelihood), 0.3846289, .01)
omxCheckCloseEnough(fivenum(m2@expectation@em.expected),
                    c(9e-05, 9.86002, 109.04181, 276.1887, 503.42342), .01)
omxCheckCloseEnough(sum(m2@expectation@em.expected), 90449.05, .01)

m2 <- mxModel(m2,
	      mxComputeIterate(steps=list(
				 mxComputeOnce('expectation', context='EM'),
				 mxComputeOnce('fitfunction', gradient=TRUE, hessian=TRUE)
				 )))
m2 <- mxRun(m2)
omxCheckCloseEnough(m2@fitfunction@result, 3221.826, .01)
omxCheckCloseEnough(fivenum(m2@output$gradient), c(-128.034, -8.294, 10.7, 25.814, 107.966), .01)
omxCheckCloseEnough(fivenum(m2@output$hessian[m2@output$hessian != 0]),
                    c(6.559, 6.559, 32.976, 83.554, 107.714), .01)

m2 <- mxModel(m2,
              mxData(observed=data, type="raw"),  # got sorted, add it again unsorted
              mxExpectationBA81(
                ItemSpec=items, EItemParam="EItemParam",ItemParam="itemParam",
                mean="mean", cov="cov",
                qpoints=31,
                scores="full"),
	      mxComputeIterate(steps=list(
				 mxComputeOnce('expectation', context='EM'),
				 mxComputeNewtonRaphson(free.set='itemParam'),
				 mxComputeOnce('expectation'),
				 mxComputeOnce('fitfunction', free.set=c("mean","cov"))
				 )))

	m2 <- mxOption(m2, "Analytic Gradients", 'Yes')
	m2 <- mxOption(m2, "Verify level", '-1')
m2 <- mxOption(m2, "Function precision", '1.0E-5')
m2 <- mxRun(m2)

#print(m2@matrices$itemParam@values)
#print(correct.mat)
omxCheckCloseEnough(m2@fitfunction@result, 6216.272, .01)
got <- cor(c(m2@matrices$itemParam@values[1:2,]),
           c(correct.mat[1:2,]))
omxCheckCloseEnough(got, .988, .01)
scores <- m2@expectation@scores.out
omxCheckCloseEnough(scores[1:5,1], c(0.6783773, 0.2848123, -0.3438632, -0.1026575, -1.0820213), .001)
omxCheckCloseEnough(scores[1:5,2], c(0.6769653, 0.6667262, 0.6629124, 0.6624804, 0.6796952), 1e-5)
omxCheckCloseEnough(scores[,1], as.vector(ability), 3.5*max(scores[,2]))
omxCheckCloseEnough(cor(c(scores[,1]), ability), .737, .01)
