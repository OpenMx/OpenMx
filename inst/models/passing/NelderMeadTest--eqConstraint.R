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
#Compare Nelder-Mead to the GD optimizer best at handling MxConstraints:
if(mxOption(NULL,"Default optimizer")!="SLSQP"){stop("SKIP")}

#The naive "soft" method only seems to work OK if EVERY vertex of the initial simplex is feasible...:
ism <- matrix(0.2,4,4) + diag(0.2,4)
colnames(ism) <- c("pred","pyellow","pgreen","pblue")
foo <- mxComputeNelderMead(iniSimplexMat=ism, nudgeZeroStarts=FALSE, xTolProx=1e-12, fTolProx=1e-8)
#foo$verbose <- 5L
plan <- omxDefaultComputePlan()
plan$steps <- list(foo,plan$steps$RE)

#Run with SLSQP:
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
	mxConstraint(Pred + Pyellow + Pgreen + Pblue - 1 == 0,name="indentifying")
)
m1run <- mxRun(m1)
summary(m1run)

m2 <- mxModel(
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
	mxConstraint(Pred + Pyellow + Pgreen + Pblue - 1 == 0,name="indentifying")
)
m2run <- mxRun(m2)
summary(m2run)

omxCheckCloseEnough(m1run$output$estimate, m2run$output$estimate, 0.01)

#Run with Nelder-Mead, with different arguments:
ism3 <- ism
ism3[4,4] <- 0.46
#^^^Note that now, it is not the case that all of the initial vertices are feasible
colnames(ism3) <- c("pred","pyellow","pgreen","pblue")
foo3 <- mxComputeNelderMead(
	iniSimplexMat=ism3, nudgeZeroStarts=FALSE, xTolProx=1e-12, fTolProx=1e-8, eqConstraintMthd="backtrack")
#foo3$verbose <- 5L
plan3 <- omxDefaultComputePlan()
plan3$steps <- list(foo3,plan3$steps$RE)
m3 <- mxModel(
	"MultinomialWithLinearConstraints",
	plan3,
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pred",name="Pred",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pyellow",name="Pyellow",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pgreen",name="Pgreen",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pblue",name="Pblue",lbound=0,ubound=1),
	mxAlgebra( -2*(43*log(Pred) + 22*log(Pyellow) + 20*log(Pgreen) + 15*log(Pblue)), name="fitfunc"),
	mxAlgebra( cbind(-2*43/Pred,-2*22/Pyellow,-2*20/Pgreen,-2*15/Pblue), name="objgrad",
						 dimnames=list(NULL,c("pred","pyellow","pgreen","pblue"))),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad",numObs=100),
	mxCI(c("pred","pyellow","pgreen","pblue")),
	mxConstraint(Pred + Pyellow + Pgreen + Pblue - 1 == 0,name="indentifying")
)
m3run <- mxRun(m3)
#The backtrack method isn't that helpful in this case:
summary(m3run)

#l1p:
foo4 <- mxComputeNelderMead(
	iniSimplexMat=ism3, nudgeZeroStarts=FALSE, xTolProx=1e-12, fTolProx=1e-8, eqConstraintMthd="l1p")
#foo3$verbose <- 5L
plan4 <- omxDefaultComputePlan()
plan4$steps <- list(foo4,plan4$steps$RE)
m4 <- mxModel(
	"MultinomialWithLinearConstraints",
	plan4,
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pred",name="Pred",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pyellow",name="Pyellow",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pgreen",name="Pgreen",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pblue",name="Pblue",lbound=0,ubound=1),
	mxAlgebra( -2*(43*log(Pred) + 22*log(Pyellow) + 20*log(Pgreen) + 15*log(Pblue)), name="fitfunc"),
	mxAlgebra( cbind(-2*43/Pred,-2*22/Pyellow,-2*20/Pgreen,-2*15/Pblue), name="objgrad",
						 dimnames=list(NULL,c("pred","pyellow","pgreen","pblue"))),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad",numObs=100),
	mxCI(c("pred","pyellow","pgreen","pblue")),
	mxConstraint(Pred + Pyellow + Pgreen + Pblue - 1 == 0,name="indentifying")
)
m4run <- mxRun(m4)
#The l1p isn't that helpful in this case, either:
summary(m4run)
#Penalized fit should be slightly greater than raw fit:
omxCheckTrue(m4run$compute$steps[[1]]$output$penalizedFit > m4run$output$fit)

#GDsearch:
foo5 <- mxComputeNelderMead(
	iniSimplexMat=ism3, nudgeZeroStarts=FALSE, xTolProx=1e-12, fTolProx=1e-8, eqConstraintMthd="GDsearch")
#foo3$verbose <- 5L
plan5 <- omxDefaultComputePlan()
plan5$steps <- list(foo5,plan5$steps$RE)
m5 <- mxModel(
	"MultinomialWithLinearConstraints",
	plan5,
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pred",name="Pred",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pyellow",name="Pyellow",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pgreen",name="Pgreen",lbound=0,ubound=1),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0.25,labels="pblue",name="Pblue",lbound=0,ubound=1),
	mxAlgebra( -2*(43*log(Pred) + 22*log(Pyellow) + 20*log(Pgreen) + 15*log(Pblue)), name="fitfunc"),
	mxAlgebra( cbind(-2*43/Pred,-2*22/Pyellow,-2*20/Pgreen,-2*15/Pblue), name="objgrad",
						 dimnames=list(NULL,c("pred","pyellow","pgreen","pblue"))),
	mxFitFunctionAlgebra(algebra="fitfunc",gradient="objgrad",numObs=100),
	mxCI(c("pred","pyellow","pgreen","pblue")),
	mxConstraint(Pred + Pyellow + Pgreen + Pblue - 1 == 0,name="indentifying")
)
m5run <- mxRun(m5)
#The GDsearch method actually gets the answer:
summary(m5run)
omxCheckCloseEnough(m1run$output$estimate, m5run$output$estimate, 1e-4)
