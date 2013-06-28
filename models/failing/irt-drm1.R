library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 10
i1 <- rpf.drm(numChoices=4)
items <- vector("list", numItems)
correct <- vector("list", numItems)
for (ix in 1:numItems) {
  items[[ix]] <- i1
  correct[[ix]] <- rpf.rparam(i1)
  correct[[ix]][3] <- 0
}
correct.mat <- simplify2array(correct)

ability <- rnorm(500)
data <- rpf.sample(ability, items, correct)

spec <- mxMatrix(name="ItemSpec", nrow=6, ncol=numItems,
         values=sapply(items, function(m) slot(m,'spec')),
         free=FALSE, byrow=TRUE)

ip.mat <- mxMatrix(name="itemParam", nrow=3, ncol=numItems,
                   values=c(1,0,0),
                   free=c(TRUE, TRUE, FALSE),
		   lbound=c(1e-6, -1e6, 0))

m2 <- mxModel(model="drm1", ip.mat, spec,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec="ItemSpec",
                ItemParam="itemParam",
		GHpoints=30),
              mxFitFunctionBA81())

m2 <- mxOption(m2, "Analytic Gradients", 'no')
if (1) {
	m2 <- mxOption(m2, "Analytic Gradients", 'yes')
	m2 <- mxOption(m2, "Verify level", '-1')
}
m2 <- mxOption(m2, "Function precision", '1.0E-5')
m2 <- mxOption(m2, "Calculate Hessian", "No")
m2 <- mxOption(m2, "Standard Errors", "No")
m2 <- mxRun(m2)

#print(m2@matrices$itemParam@values)
#print(correct.mat)
got <- cor(c(m2@matrices$itemParam@values),
           c(correct.mat))
omxCheckCloseEnough(got, .988, .01)
ability <- scale(ability)
omxCheckCloseEnough(m2@output$ability[,1], as.vector(ability), 3.5*max(m2@output$ability[,2]))
omxCheckCloseEnough(cor(c(m2@output$ability[,1]), ability), .80, .01)
q()
