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


require(OpenMx)

set.seed(99)
data <- data.frame(x=rnorm(10), y=as.numeric(sample(1:2, size=10, replace=TRUE)))

for (toSortData in c(TRUE,FALSE)) {
  model <- mxModel('Test Definition Row Sorting',
                   mxMatrix('Full', 1, 2, values=0, free=TRUE, name='TwoMeans'),
                   mxMatrix('Full', 1, 1, values=1, free=FALSE, name='Var'),
                   mxAlgebra(TwoMeans[1, data.y], 'expMean'),
                   mxData(data, 'raw', sort=toSortData),
                   mxFitFunctionML(),
                   mxExpectationNormal('Var', 'expMean', dimnames='x')
  )
  
  fit <- mxRun(model)
  
  fit.means <- omxGetParameters(fit)
  obs.means <- c(mean(data[data$y==1,'x']), mean(data[data$y==2,'x']))
  
  omxCheckCloseEnough(fit.means, obs.means, 0.001)
  
  perRow <- sapply(1:10, function(row) mxEval(expMean, fit, defvar.row=row, compute=T))
  group <- apply(cbind(abs(perRow - fit.means[1]),
                       abs(perRow - fit.means[2])), 1, order)[1,]
  omxCheckEquals(data$y, group)
}
