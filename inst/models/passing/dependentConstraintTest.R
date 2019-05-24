#
#   Copyright 2007-2018 by the individuals mentioned in the source code history
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

#Currently, only NPSOL can deal with non-identical but linearly dependent equality constraints.

library(OpenMx)
if(mxOption(NULL,"Default optimizer")=="CSOLNP"){stop("SKIP")}

if(mxOption(NULL,"Default optimizer")=="NPSOL"){
	m1 <- mxModel(
		"MultinomialWithLinearConstraints",
		mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pred",name="Pred",lbound=0,ubound=1),
		mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pyellow",name="Pyellow",lbound=0,ubound=1),
		mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pgreen",name="Pgreen",lbound=0,ubound=1),
		mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pblue",name="Pblue",lbound=0,ubound=1),
		mxAlgebra( -2*(43*log(Pred) + 22*log(Pyellow) + 20*log(Pgreen) + 15*log(Pblue)), name="fitfunc"),
		mxAlgebra( cbind(-2*43/Pred,-2*22/Pyellow,-2*20/Pgreen,-2*15/Pblue), name="objgrad",
							 dimnames=list(NULL,c("pred","pyellow","pgreen","pblue"))),
		mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad",numObs=100),
		mxCI(c("pred","pyellow","pgreen","pblue")),
		mxConstraint(Pred + Pyellow + Pgreen + Pblue - 1 == 0,name="indentifying"),
		mxConstraint(2*(Pred + Pyellow + Pgreen + Pblue) - 2 == 0,name="indentifying2")
	)
	m1run <- mxRun(m1)
	summary(m1run)
	omxCheckCloseEnough(coef(m1run),c(0.43,0.22,0.2,0.15),1e-7)
}

if(mxOption(NULL,"Default optimizer")=="SLSQP"){
	plan <- omxDefaultComputePlan()
	plan$steps <- list(GD=plan$steps$GD)
	plan$steps$GD <- mxComputeNelderMead(nudgeZeroStarts=FALSE, eqConstraintMthd="GDsearch")
	
	m1 <- mxModel(
		"MultinomialWithLinearConstraints",
		plan,
		mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pred",name="Pred",lbound=0,ubound=1),
		mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pyellow",name="Pyellow",lbound=0,ubound=1),
		mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pgreen",name="Pgreen",lbound=0,ubound=1),
		mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pblue",name="Pblue",lbound=0,ubound=1),
		mxAlgebra( -2*(43*log(Pred) + 22*log(Pyellow) + 20*log(Pgreen) + 15*log(Pblue)), name="fitfunc"),
		mxAlgebra( cbind(-2*43/Pred,-2*22/Pyellow,-2*20/Pgreen,-2*15/Pblue), name="objgrad",
							 dimnames=list(NULL,c("pred","pyellow","pgreen","pblue"))),
		mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad",numObs=100),
		mxCI(c("pred","pyellow","pgreen","pblue")),
		mxMatrix(type="Unit",nrow=1,ncol=4,dimnames=list(NULL,c("pred","pyellow","pgreen","pblue")),name="jac"),
		mxAlgebra(2%x%jac, name="jac2", dimnames=list(NULL,c("pred","pyellow","pgreen","pblue"))),
		mxConstraint(Pred + Pyellow + Pgreen + Pblue - 1 == 0,name="indentifying",jac="jac"),
		mxConstraint(2*(Pred + Pyellow + Pgreen + Pblue) - 2 == 0,name="indentifying2",jac="jac2")
	)
	omxCheckWarning(
		mxRun(m1),
		"counted 2 equality constraints, but equality-constraint Jacobian is apparently rank 1 at the start values; Nelder-Mead will not work correctly unless equality constraints are linearly independent (this warning may be spurious if there are non-smooth equality constraints)"
	)
	
}
