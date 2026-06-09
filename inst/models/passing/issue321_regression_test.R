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
