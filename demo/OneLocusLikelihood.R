#
#   Copyright 2007-2012 The OpenMx Project
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

# -----------------------------------------------------------------------------
# Program: OneLocusLikelihood.R  
# Author: Michael Neale
# Date: 2009.11.23 
#
# ModelType: Locus Likelihood
# DataType: Frequencies
# Field: Human Behavior Genetics
#
# Purpose: 
#      One Locus Likelihood model to estimate allele frequencies
#      Bernstein data on ABO blood-groups
#      c.f. Edwards, AWF (1972)  Likelihood.  Cambridge Univ Press, pp. 39-41
#
# RevisionHistory:
#      Hermine Maes -- 2010.02.22 updated & reformatted
#      Ross Gore -- 2011.06.15 added Model, Data & Field
# -----------------------------------------------------------------------------

require(OpenMx)
# Load Library
# -----------------------------------------------------------------------------
    
OneLocusModel <- mxModel("OneLocus",
    mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=c(.3333), name="P"),
    mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=c(.3333), name="Q"),
    mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=c(.3333), name="R"),
    # Matrices for allele frequencies, p, q and r
    # -------------------------------------
    mxConstraint(P+Q+R == 1, name = "EqualityConstraint"),
    mxMatrix( type="Full", nrow=4, ncol=1, values=c(211,104,39,148), name="ObservedFreqs"),
	# Matrix of observed data
    # -------------------------------------
    mxAlgebra( expression=rbind(P*(P+2*R), Q*(Q+2*R), 2*P*Q, R*R), name="ExpectedFreqs"),
	# Algebra for predicted proportions
    # -------------------------------------
    mxAlgebra( expression=-(sum(log(ExpectedFreqs) * ObservedFreqs)), name="NegativeLogLikelihood"),
    # Algebra for -logLikelihood
    # -------------------------------------
    mxAlgebraObjective("NegativeLogLikelihood") 
    # User-defined objective
     # -------------------------------------
)
# Create an MxModel object
# -----------------------------------------------------------------------------

OneLocusFit <- mxRun(OneLocusModel)
OneLocusFit@matrices
OneLocusFit@algebras


estimates <- c(
    OneLocusFit@matrices$P@values, 
    OneLocusFit@matrices$Q@values, 
    OneLocusFit@matrices$R@values, 
    OneLocusFit@algebras$NegativeLogLikelihood@result)
Mx1Estimates<-c(0.2931,0.1552,0.5517,627.851)
omxCheckCloseEnough(estimates,Mx1Estimates,.01)
# Compare OpenMx estimates to original Mx's estimates
# -----------------------------------------------------------------------------

