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

#Two locus Likelihood example
#Author: Mike Neale
#Date: Nov 23 2009

#     Bernstein data on ABO blood-groups
#     c.f. Edwards, AWF (1972)  Likelihood.  Cambridge Univ Press, pp. 39-41
#
twolocus<-mxModel("twolocus", mxFitFunctionAlgebra("NegativeLogLikelihood"), 
				
				mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=c(.3333), name="P"),
				mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=c(.3333), name="S"),
				mxMatrix("Full", nrow=4, ncol=1, values=c(212,103,39,148),name="ObservedFreqs"),

				mxAlgebra(1-P, name="Q"),
				mxAlgebra(1-S, name="T"),
				mxAlgebra(rbind ((P*P+2*P*Q)*T*T, (Q*Q)*(S*S+2*S*T), (P*P+2*P*Q)*(S*S+2*S*T), (Q*Q)*(T*T)), name="ExpectedFreqs"),
				mxAlgebra(-(sum( log(ExpectedFreqs) * ObservedFreqs )), name = "NegativeLogLikelihood")
			)

#run the model
run<-mxRun(twolocus, suppressWarnings=TRUE)
run$matrices
run$algebras


# Compare OpenMx estimates to original Mx's estimates
estimates<-c(run$matrices$P$values,run$matrices$S$values,run$algebras$NegativeLogLikelihood$result)
Mx1Estimates<-c(0.2929,0.1532,646.972)
omxCheckCloseEnough(estimates,Mx1Estimates,.01)

