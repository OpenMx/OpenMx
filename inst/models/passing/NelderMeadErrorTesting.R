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

#This test script doesn't use any of the GD optimizers, so there's no reason to run it when all 3 are the default:
if(mxOption(NULL,"Default optimizer")!="CSOLNP"){stop("SKIP")}


omxCheckError(
	mxComputeNelderMead(unrecognizedArgument=3),
	"mxComputeNelderMead() does not accept values for the '...' argument"
)

omxCheckError(
	mxComputeNelderMead(nudgeZeroStarts="foo"),
	"unrecognized character string provided as argument 'nudgeZeroStarts'"
)

omxCheckError(
	mxComputeNelderMead(alpha=-1),
	"reflection coefficient 'alpha' must be positive"
)

omxCheckError(
	mxComputeNelderMead(betao=-1),
	"contraction coefficients 'betao' and 'betai' must both be within unit interval (0,1)"
)
omxCheckError(
	mxComputeNelderMead(betao=2),
	"contraction coefficients 'betao' and 'betai' must both be within unit interval (0,1)"
)
omxCheckError(
	mxComputeNelderMead(betai=-1),
	"contraction coefficients 'betao' and 'betai' must both be within unit interval (0,1)"
)
omxCheckError(
	mxComputeNelderMead(betai=2),
	"contraction coefficients 'betao' and 'betai' must both be within unit interval (0,1)"
)

omxCheckError(
	mxComputeNelderMead(alpha=1,gamma=0.9),
	"if positive, expansion coefficient 'gamma' must be greater than reflection coefficient 'alpha'"
)

omxCheckError(
	mxComputeNelderMead(sigma=1),
	"shrink coefficient 'sigma' must be less than 1.0"
)

omxCheckError(
	mxComputeNelderMead(iniSimplexType="foo"),
	"'foo' should be one of 'regular', 'right', 'smartRight', and 'random'"
)

omxCheckError(
	mxComputeNelderMead(degenLimit=-1),
	"'degenLimit' must be within interval [0,pi]"
)
omxCheckError(
	mxComputeNelderMead(degenLimit=-1),
	"'degenLimit' must be within interval [0,pi]"
)

omxCheckError(
	mxComputeNelderMead(stagnCtrl=10),
	"'stagnCtrl' must be an integer vector of length 2"
)

omxCheckError(
	mxComputeNelderMead(ineqConstraintMthd="foo"),
	"'foo' should be one of 'soft' and 'eqMthd'"
)
omxCheckError(
	mxComputeNelderMead(eqConstraintMthd="foo"),
	"'foo' should be one of 'GDsearch', 'soft', 'backtrack', and 'l1p'"
)

##################################################

set.seed(1611150)
x <- matrix(rnorm(1000,sd=2))
colnames(x) <- "x"
m <- mxModel(
	"mod",
	mxData(observed=x,type="raw"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=0,labels="mu",name="Mu"),
	mxMatrix(type="Full",nrow=1,ncol=1,free=T,values=4,labels="sigma2",name="Sigma2",lbound=0.0001),
	mxExpectationNormal(covariance="Sigma2",means="Mu",dimnames=c("x")),
	mxAlgebra(sqrt(Sigma2),name="Sigma"),
	mxFitFunctionML()
)
plan <- omxDefaultComputePlan()
plan$steps$GD <- mxComputeNelderMead()
plan2 <- plan

plan$steps$GD$alpha <- -1
omxCheckError(mxRun(mxModel(m,plan)), "reflection coefficient 'alpha' must be positive")
plan <- plan2

plan$steps$GD$betao <- -1
omxCheckError(mxRun(mxModel(m,plan)), "contraction coefficients 'betao' and 'betai' must both be within unit interval (0,1)")
plan <- plan2

plan$steps$GD$betai <- 1
omxCheckError(mxRun(mxModel(m,plan)), "contraction coefficients 'betao' and 'betai' must both be within unit interval (0,1)")
plan <- plan2

plan$steps$GD$gamma <- 0.9
omxCheckError(mxRun(mxModel(m,plan)), "if positive, expansion coefficient 'gamma' must be greater than reflection coefficient 'alpha'")
plan <- plan2

plan$steps$GD$iniSimplexType <- "foo"
omxCheckError(mxRun(mxModel(m,plan)), "unrecognized character string provided for Nelder-Mead 'iniSimplexType'")
plan <- plan2

plan$steps$GD$degenLimit <- -1
omxCheckError(mxRun(mxModel(m,plan)), "'degenLimit' must ge within interval [0,pi]")
plan <- plan2

plan$steps$GD$stagnCtrl <- 10L
omxCheckError(mxRun(mxModel(m,plan)), "'stagnCtrl' must be an integer vector of length 2")
plan <- plan2

plan$steps$GD$xTolProx <- -1
plan$steps$GD$fTolProx <- -1
omxCheckWarning(mxRun(mxModel(m,plan)), "both 'xTolProx' and 'fTolProx' are non-positive; 'fTolProx' will be assigned a value of 1e-14")
plan <- plan2

plan$steps$GD$ineqConstraintMthd <- "foo"
omxCheckError(mxRun(mxModel(m,plan)), "unrecognized character string provided for Nelder-Mead 'ineqConstraintMthd'")
plan <- plan2

plan$steps$GD$eqConstraintMthd <- "foo"
omxCheckError(mxRun(mxModel(m,plan)), "unrecognized character string provided for Nelder-Mead 'eqConstraintMthd'")
plan <- plan2

ism <- matrix(0,1,1)
plan$steps$GD$iniSimplexMat <- ism
omxCheckError(mxRun(mxModel(m,plan)), "'iniSimplexMat' has 1 columns, but 2 columns expected")
plan <- plan2

ism <- matrix(c(0,4,-1,4,0,5,1,1),4,2)
colnames(ism) <- c("mu","sigma2")
plan$steps$GD$iniSimplexMat <- ism
plan$steps$GD$fTolProx <- 1e-8
plan$steps$GD$xTolProx <- 1e-8
omxCheckWarning(mxRun(mxModel(m,plan)), "'iniSimplexMat' has 4 rows, but 3 rows expected; extraneous rows will be ignored")
plan <- plan2

ism <- matrix(c(0,4,-1,4,0,5),3,2,byrow=T)
colnames(ism) <- c("um","sigma2")
plan$steps$GD$iniSimplexMat <- ism
omxCheckError(mxRun(mxModel(m,plan)), "error in mapping column names of 'iniSimplexMat' to free-parameter labels")
plan <- plan2

ism <- matrix(c(0,4,-1,5),2,2,byrow=T)
colnames(ism) <- c("mu","sigma2")
plan$steps$GD$iniSimplexMat <- ism
plan$steps$GD$fTolProx <- 1e-8
plan$steps$GD$xTolProx <- 1e-8
omxCheckWarning(mxRun(mxModel(m,plan)), "'iniSimplexMat' has 2 rows, but 3 rows expected; omitted rows will be generated randomly")
plan <- plan2

ism <- matrix(c(2,4,1,4,3,5),3,2,byrow=T)
colnames(ism) <- c("mu","sigma2")
plan$steps$GD$iniSimplexMat <- ism
#plan$steps$GD$verbose <- 5L
omxCheckError(
	mxRun(mxModel(m,plan,mxConstraint(Mu<0))),
	"The job for model 'mod' exited abnormally with the error message: initial simplex is not feasible; specify it differently, try different start values, or use mxTryHard()"
)
plan <- plan2

ism <- matrix(c(-2,4,-1,4,-3,5),3,2,byrow=T)
colnames(ism) <- c("mu","sigma2")
plan$steps$GD$iniSimplexMat <- ism
#plan$steps$GD$verbose <- 5L
omxCheckError(
	mxRun(mxModel(m,plan,mxConstraint(Mu>0))),
	"The job for model 'mod' exited abnormally with the error message: initial simplex is not feasible; specify it differently, try different start values, or use mxTryHard()"
)
plan <- plan2

ism <- matrix(c(0,4,-1,4,0,5),3,2,byrow=T)
plan$steps$GD$iniSimplexMat <- ism
omxCheckError(mxRun(mxModel(m,plan)), "'iniSimplexMat' has 0 column names, but 2 column names expected")
plan <- plan2
