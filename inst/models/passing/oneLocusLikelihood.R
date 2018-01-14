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

#Single locus Likelihood example
#Author: Mike Neale
#Date: Nov 23 2009

#     Bernstein data on ABO blood-groups
#     c.f. Edwards, AWF (1972)  Likelihood.  Cambridge Univ Press, pp. 39-41
#
onelocus<-mxModel("onelocus", mxFitFunctionAlgebra("NegativeLogLikelihood"), 
				
				mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=c(.3333), name="P"),  # P, freq of allele 1
				mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=c(.3333), name="Q"),  # Q, freq of allele 2
				mxMatrix("Full", nrow=1, ncol=1, free=TRUE, values=c(.3333), name="R"),  # R, freq of allele 3
				mxMatrix("Unit", nrow=1, ncol=1, name="I"),                              # 1.0 constant for equality constraint
				mxMatrix("Full", nrow=4, ncol=1, values=c(212,103,39,148),name="ObservedFreqs"), # Data

				mxAlgebra(rbind(P*(P+2*R), Q*(Q+2*R), 2*P*Q, R*R), name="ExpectedFreqs"), # Predicted proportions
				mxAlgebra(-(sum( log(ExpectedFreqs) * ObservedFreqs )), name = "NegativeLogLikelihood"),  # here is -log Likelihood
				
				mxConstraint(P + Q + R == 1, "equalityConstraint")   # This is the equality constraint
			)

#run the model
run<-mxRun(onelocus)
#run$matrices
#run$algebras

omxCheckCloseEnough(run$algebras$NegativeLogLikelihood$result, 627.028, .1)

# Compare OpenMx estimates to original Mx's estimates
estimates<-c(run$matrices$P$values,run$matrices$Q$values,run$matrices$R$values)
Mx1Estimates<-c(0.2945,0.1540,0.5515)
omxCheckCloseEnough(estimates,Mx1Estimates,.01)

if (mxOption(NULL, 'Default optimizer') != "CSOLNP") {
	onelocus <- mxModel(onelocus,
		    mxConstraint(P + Q + R == 1, "redundent"))
	run <- mxRun(onelocus)
	omxCheckCloseEnough(run$algebras$NegativeLogLikelihood$result, 627.028, .1)
	estimates<-c(run$matrices$P$values,run$matrices$Q$values,run$matrices$R$values)
	omxCheckCloseEnough(estimates,Mx1Estimates,.01)
}
