#
#   Copyright 2007-2026 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

library(OpenMx)

# Skip test if the default optimizer is not SLSQP
if(mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")

set.seed(12345)
df = MASS::mvrnorm(n=100, mu=c(0,0), Sigma=matrix(c(1,.5,.5,1), ncol=2, nrow=2))
colnames(df) <- c("x", "y")
df = data.frame(df)

# Fit a simple regression model on covariance data
model_cov = mxModel("Simple Regression Covariance Data", 
    type = "RAM",
    manifestVars = c("x", "y"),
    mxPath(from = "x", to = "y", arrows = 1, free = TRUE, values = 0.5, labels = "beta1"),
    mxPath(from = c("x", "y"), arrows = 2, free = TRUE, values = 1, labels = c("varx", "sigma2")),
    mxPath(from = "one", to = c("x", "y"), arrows = 1, free = TRUE, values = 0, labels = c("meanx", "beta0")),
    mxData(observed=cov(df), type="cov", means=apply(df, 2, mean), numObs=nrow(df))
)
fit_cov = mxRun(model_cov)
sum_cov = summary(fit_cov)

# Since the model is fully saturated (df = 0), the chi-square fit statistic
# should be exactly 0 (within floating-point precision) under both the default 
# summary formula and the run-based reference models.
omxCheckCloseEnough(sum_cov$Chi, 0.0, 1e-4)

# Fit with explicit reference models and verify they match
ref_cov = mxRefModels(fit_cov, run=TRUE)
sum_cov_ref = summary(fit_cov, refModels = ref_cov)
omxCheckCloseEnough(sum_cov_ref$Chi, 0.0, 1e-4)

# Confirm they are equal to each other
omxCheckCloseEnough(sum_cov$Chi, sum_cov_ref$Chi, 1e-6)


# More tests, based on Mike Cheung's syntax: #####

## mxModel using the parameter estimates as starting values
matrA <- mxMatrix(
	type="Full", nrow=2, ncol=2,
	free=c(F,F,T,F), values=0,
	labels=c(NA,NA,"beta1",NA), byrow=TRUE, name="A")
matrS <- mxMatrix(
	type="Symm", nrow=2, ncol=2,
	free=c(T,F,F,T), values=0,
	labels=c("varx",NA,NA,"sigma2"), byrow=TRUE, name="S")
matrF <- mxMatrix(type="Iden", nrow=2, ncol=2, name="F")
matrM <- mxMatrix(
	type="Full", nrow=1, ncol=2,
	free=c(T,T), values=0,
	labels=c("meanx","beta0"), name="M")
expRAM <- mxExpectationRAM("A","S","F","M", dimnames=c("x","y"))
funML <- mxFitFunctionML()
## Raw data as inputs
uniModel1 <- mxModel(
	"Simple Regression with raw data",
	mxData(observed=df, type="raw"),
	matrA, matrS, matrF, matrM, expRAM, funML)
fit1 <- mxRun(uniModel1)
( smm1 <- summary(fit1, refModels=mxRefModels(fit1, run=TRUE)) )

## Covariance matrix as inputs
uniModel2 <- mxModel(
	"Simple Regression with covariance matrix",
	mxData(
		observed=cov(df), type="cov",
		means=apply(df, 2, mean),
		numObs=nrow(df)),
	matrA, matrS, matrF, matrM, expRAM, funML)
fit2 <- mxRun(uniModel2)
( smm2 <- summary(fit2) )

# Chi-square statistics should match very closely:
omxCheckCloseEnough(smm1$Chi,smm2$Chi,1e-12)
