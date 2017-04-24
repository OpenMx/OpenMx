# ===========
# = HISTORY =
# ===========
# 2017-04-14 04:57PM TBATES
# 	update for   mxExpectationRAM("A", "S", "F", "M", thresholds="thresholds"),mxFitFunctionML()
# 	fix "Function precision" bug
# has a bug (no dimnames) at 104
library(OpenMx)
library(mvtnorm)

# define sample size
n <- 250
# define loadings for simulation model when covariate=0
loadings1 <- c(1, 1, 1, 1, 1)
# define change in loadings dependent on covariate
loadings2 <- c(0, 0, 0, 0, -.5)
# define residual variance terms
resid <- c(.4, .4, .4, .4, .4)
# define definition variable
g <- runif(n, 0, 2)
# simulate data for testing
data <- matrix(NA, n, 6)
for (i in 1:n){
	l <- loadings1 + loadings2 * g[i]
	e <- matrix(0, 5, 5)
	diag(e) <- resid
	data[i,1:5] <- rnorm(1)*l + rmvnorm(1,rep(0,5),e)
}
data <- data.frame(data)
names(data) <- c(paste("x",1:5, sep=""),"g")
# constrain the data to categories
data[data < 0] <- 0
data[data > 0] <- 1
# add the definition variable
data[, 6] <- g

# view the correlation matrix of the data
cor(data)

# ================
# = Begin OpenMx =
# ================

## Let's try a threshold model without a definition variable##
A <- mxMatrix("Full", 6, 6, name="A")
A$values[1:5, 6] <- 1
A$free[1:5, 6]   <- TRUE
A$labels[1:5, 6] <- paste("l",1:5,sep="")

S <- mxMatrix("Symm", 6, 6, name="S")
diag(S$values) <- 1
diag(S$labels) <- c(paste("e",1:5,sep=""), "var.F")

F <- mxMatrix("Full", 5, 6, dimnames=list(names(data)[1:5], c(names(data)[1:5], "F")), name="F")
diag(F$values) <- 1

M <- mxMatrix("Zero", 1, 6, name="M")

# columns are variables; there are four thresholds for each variable
T <- mxMatrix("Full", 1, 5, TRUE, 0, dimnames = list(c(), names(data)[1:5]), name="thresholds")

catModel <- mxModel("Categorical Test", A, S, F, M, T,
	mxData(data, type="raw"),
  mxExpectationRAM("A", "S", "F", "M", thresholds="thresholds"),
	mxFitFunctionML()
)

catResults <- mxRun(catModel)

summary(catResults)

# item difficulties
catResults$output$estimate[6:10]/catResults$output$estimate[1:5]

# item discriminations
catResults$output$estimate[1:5]


## Let's try a threshold model with a definition variable##

# here's the trick. i created a dummy factor called "hack", 
# and labeled it's path on x5 with the name of the definition variable
# there's a flashier way involving mxAlgebra statements, but the required substitution isn't working
# look for "square bracket substitution" in the future
CA <- mxMatrix("Full", 7, 7, name="CA")
CA$values[c(1:5,7), 6] <- 1
CA$free[c(1:5,7), 6]   <- TRUE
CA$labels[c(1:5,7), 6] <- c(paste("l",1:5,sep=""), "dif")

CA$labels[5, 7] <-"data.g"
#trick <- mxMatrix("Full", 1, 2, TRUE, c(1,0), labels=c("alpha","beta"))
#def   <- mxAlgebra(alpha + beta*data.g, name = "l5")

CS <- mxMatrix("Symm", 7, 7, name="CS")
diag(CS$values) <- c(rep(1,6),0)
diag(CS$labels) <- c(paste("e",1:5,sep=""), "var.F", "var.hack")

CF <- mxMatrix("Full", 5, 7, dimnames=list(names(data)[1:5], c(names(data)[1:5], "F", "hack")), name="CF")
diag(CF$values) <- 1

# =================================================
# = FAILS to RUN DUE TO LACKING DIMNAMES ON MEANS =
# =================================================
CM <- mxMatrix("Zero", 1, 7, name = "CM",)

# columns are variables; there are four thresholds for each variable
CT <- mxMatrix("Full", 1, 5, TRUE, 0, 
    labels = paste("t", 1:5, sep=""), 
    dimnames = list(c(), names(data)[1:5]), 
    name="thresholds")

catModel2 <- mxModel("Categorical Test with Def", CA, CS, CF, CM, CT, 
    mxData(data, type="raw"),
	  mxExpectationRAM("CA", "CS", "CF", "CM", thresholds="thresholds"),
		mxFitFunctionML()		
)

catModel2 <- mxOption(catModel2, "Function precision", 1e-9)

catResults2 <- mxRun(catModel2)

summary(catResults2)


# item difficulties
catResults2$output$estimate[6:10]/catResults2$output$estimate[1:5]

# item discriminations
catResults2$output$estimate[1:5]

catResults2$output
