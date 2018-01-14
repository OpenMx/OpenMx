#
#   Copyright 2007-2018 The OpenMx Project
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

#   A very big thank you to Tim Bates for writing this script.

require(OpenMx)
#options(error = utils::recover)
#options(warn = 2)

stdErr <- function(x) {
	sqrt(var(x)/length(x))
}

set.seed(200); 
known <- rnorm(1000,mean=0,sd=1)
known <- as.data.frame(known); 
names(known)<-c("test")
stdErr(known$test) # 0.03145627

manifests <- names(known)
latents <- c()
factorModel <- mxModel("test SE", type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=manifests, arrows=2, values=.1), # variance
      mxPath(from="one", to=manifests),
      mxData(known, type="raw")
)
factorModel <- mxOption(factorModel, "Standard Errors", "Fnord")
ignore <- omxCheckError(mxRun(factorModel),
			  "mxOption 'Standard Errors' must be either 'Yes' or 'No'");
factorModel <- mxOption(factorModel, "Standard Errors", "Yes")
fit <- mxRun(factorModel)
test.summary <- summary(mxRun(factorModel, suppressWarnings=TRUE))

omxCheckCloseEnough(test.summary$parameters$Std.Error[[2]], stdErr(known$test), 0.001)
