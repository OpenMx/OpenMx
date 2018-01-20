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


library(OpenMx)

dn <- paste("m",1:2, sep="")
dataCov <- matrix(c(1,.2,.2,1.5), nrow=2, dimnames=list(dn,dn))
dataMeans <- c(-.2, .3)
names(dataMeans) <- dn

n2 <- mxModel("normal2",
        mxData(observed=dataCov, type="cov",
               means=dataMeans, numObs=35),
        mxMatrix(name="cov", "Symm", 2, 2, free=T, values = c(1, 0, 1),
                 labels = c("var1", "cov12", "var2"), dimnames=list(dn,dn)),
        mxMatrix(name="mean", "Full", 1, 2, free=T, labels=dn, dimnames=list(NULL,dn)),
        mxFitFunctionML(),
        mxExpectationNormal("cov", "mean"))

plan <- mxComputeSequence(list(
  mxComputeNewtonRaphson(),
  mxComputeOnce('fitfunction', 'information', 'hessian'),
  mxComputeStandardError(),
  mxComputeHessianQuality(),
  mxComputeReportDeriv()))

for (retry in 1:2) {
  if (retry == 2) n2 <- mxModel(n2, plan)
  n2Fit <- mxRun(n2)
  
  omxCheckCloseEnough(n2Fit$output$fit, 81.216, .01)
  omxCheckCloseEnough(n2Fit$cov$values, (34/35) * dataCov, 1e-4)
  omxCheckCloseEnough(c(n2Fit$mean$values), dataMeans, 1e-4)
  omxCheckCloseEnough(log(n2Fit$output$conditionNumber), 1.63, .2)
  omxCheckCloseEnough(c(n2Fit$output$standardErrors),
                      c(0.232, 0.203, 0.348, 0.166, 0.204), .01)
}
