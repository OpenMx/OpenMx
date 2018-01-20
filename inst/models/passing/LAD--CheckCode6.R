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

# Script (by Rob K.) to check whether both NPSOL and CSOLNP return status 6 ("Code Red").  This is accomplished
# through a naive specification of a simple least-absolute-deviations (LAD) regression problem, in parametric form
# by assuming the residuals are iid LaPlace, via mxFitFunctionRow().  This can reliably produce a Code Red
# since neither optimizer will be able to zero the gradient.

require(OpenMx)
#mxOption(NULL, "Default optimizer", "CSOLNP")

#Response variable y:
y <- matrix(c(3.86735442126894,3.21609311807948,1.6681246111281,3.54171497693329,3.02206567904312,2.40194706094571,
              4.00354871075935,3.50175679256405,3.92466558500823,4.26190374144865,2.68136931922448,3.81633303744976,
              7.63868421858379,4.42267196614715,8.92858721732324,5.23749096355528,5.79233361162228,4.68998250109233,
              4.6116111993746,5.56825200215576,3.96879794025251),dimnames=list(NULL,"y"))
#Predictor variable x:
x <- matrix(c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2,2,2,2),dimnames=list(NULL,"x"))
#Data matrix:
dat <- cbind(x,y)


Laplace_rgsn_mod1 <- mxModel(
  "LaplaceReg",
  mxData(dat, type="raw"),
  mxMatrix( type="Full", nrow=1, ncol=1, 
            free=T, 
            values=5, 
            name="A", labels="a" ), #<--Intercept
  mxMatrix( type="Full", nrow=1, ncol=1, 
            free=TRUE, 
            values=0, 
            name="B", labels="b" ), #<--Slope
  mxMatrix( type="Full", nrow=1, ncol=1, 
            free=TRUE,
            values=5, 
            name="lambda", labels="lambdapar", lbound=0.0001), #<--Residual dispersion parameter
  mxMatrix(type="Full", nrow=1, ncol=1, free=F,
           name="xdef", labels="data.x"), #<--x is definition variable.
  mxAlgebra(A + B*xdef, name="yhat", dimnames=list(NULL,"y")), #<--yhat
  mxAlgebra((log(2*lambda) + (abs(filteredDataRow - yhat)/lambda) 
  ), name="rowAlgebra"), #<--Negative loglikelihood for 1 observation
  mxAlgebra(2*sum(rowResults), name="reduceAlgebra"), #<--Full-sample deviance
  mxFitFunctionRow(rowAlgebra='rowAlgebra',
                 reduceAlgebra='reduceAlgebra',
                 dimnames=c('y'))
)

Laplace_rgsn_fit1 <- mxRun(Laplace_rgsn_mod1)
summary(Laplace_rgsn_fit1)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$fit, 64.55, .02)
omxCheckCloseEnough(Laplace_rgsn_fit1$output$estimate, c(3.21, 1.01, 0.85), .03)

print(Laplace_rgsn_fit1$compute$steps[[2]]$output$gradient)
if (0) {
  # We cannot detect this problem reliably.
  omxCheckEquals(Laplace_rgsn_fit1$output$status$code,6)
  omxCheckEquals(summary(Laplace_rgsn_fit1)$npsolMessage, "The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)")
}
Laplace_rgsn_fit2 <- mxTryHard(Laplace_rgsn_fit1)
