# library(mice)
# install.packages("Amelia")
library(Amelia)

# 1. make a variable with missingness
demoOneFactor$x = demoOneFactor$x2
demoOneFactor$x[as.logical(rbinom(n=500,1,.2))] = NA

# 2. make and run a model x <- x1
m1 <- mxModel("m1", type="RAM",
	manifestVars = c("x", "x1"),
	# Factor loadings
	mxPath("x1", to = "x"),
	mxPath(c("x1", "x"), arrows = 2), # manifest residuals 
	mxPath("one", to = c("x1", "x")), # manifest means
	mxData(demoOneFactor, type = "raw")
)
m1 = umxRun(m1, setLabels = T, setValues = T)
summary(m1)$parameters

# 3. impute some data and mouse it

# imp = mice(demoOneFactor, m = 5)
imp = amelia(demoOneFactor, m = 5)

for (i in seq) {
	imp$imputations$imp1
}

m2 = with(data= imp, expr = m1)

summary(m2)$parameters
# Error: No glance method for objects of class MxRAMModel

# This error is not desirable :-)