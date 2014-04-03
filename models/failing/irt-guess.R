library(OpenMx)
library(rpf)

set.seed(9)

numItems <- 8
i1 <- rpf.drm(numChoices=2)
i2 <- rpf.drm(numChoices=2, factors=2)
items <- vector("list", numItems)
for (ix in seq(1,numItems,2)) items[[ix]] <- i1
for (ix in seq(2,numItems,2)) items[[ix]] <- i2
correct <- lapply(items, rpf.rparam)
correct.mat <- lapply(correct, function(p) if (length(p)==4) c(p,NA) else p)
correct.mat <- simplify2array(correct.mat)

numPeople <- 1000
ability <- rbind(rnorm(numPeople), rnorm(numPeople))
data <- rpf.sample(ability, items, correct)

spec <- mxMatrix(name="ItemSpec", nrow=6, ncol=numItems,
         values=sapply(items, function(m) slot(m,'spec')),
         free=FALSE, byrow=TRUE)

ip.mat <- mxMatrix(name="itemParam", nrow=5, ncol=numItems,
                   values=c(1,0,.5,1,NA,
                            1,1.4,0,.5,1),
                   free=c(rep(TRUE,3), FALSE, FALSE,
                          rep(TRUE,4), FALSE),
		   lbound=c(1e-6, -1e6, 1e-6, 0, NA,
		     1e-6, 1e-6, -1e6, 1e-6, 0),
		   ubound=c(1e3, 1e6, 1-1e-3, 1, NA,
		     1e3, 1e3, 1e6, 1-1e-3, 1))

m2 <- mxModel(model="guess", ip.mat, spec,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(
                ItemSpec="ItemSpec",
                ItemParam="itemParam",
		qpoints=19),  # use less points since accuracy is bad anyway
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

#print(m2$matrices$itemParam$values)
#print(correct.mat)
ipv <- m2$matrices$itemParam$values
got <- cor(c(ipv[!is.na(ipv)]),
           c(correct.mat[!is.na(correct.mat)]))
omxCheckCloseEnough(got, .492, .01)   # removed prior so accuracy is bad
#ability <- scale(ability)
#omxCheckCloseEnough(m2$output$ability[,1], as.vector(ability), 4*max(m2$output$ability[,2]))
#omxCheckCloseEnough(cor(c(m2$output$ability[,1]), ability), .80, .01)
