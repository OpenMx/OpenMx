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

#No need to run this test with other than the on-load default GD optimizer:
if(mxOption(NULL,"Default optimizer")!="CSOLNP"){stop("SKIP")}

#Ordinal Data test, based on poly3dz.mx (as in models/passing/OrdinalTest.R):

# Data
nthresh1 <- 1
nthresh2 <- 12	
cnames <- c("t1neur1", "t1mddd4l", "t2neur1", "t2mddd4l")
data <- suppressWarnings(try(read.table("data/mddndzf.dat", na.string=".", col.names=cnames)))
if (is(data, "try-error")) data <- read.table("../passing/data/mddndzf.dat", na.string=".", col.names=cnames)
data[,c(1,3)] <- mxFactor(data[,c(1,3)], c(0 : nthresh2))
data[,c(2,4)] <- mxFactor(data[,c(2,4)], c(0 : nthresh1))

diff <- nthresh2 - nthresh1
nvar <- 4

Mx1Threshold <- rbind(
	c(-1.9209, 0.3935, -1.9209, 0.3935),
	c(-0.5880, 0    , -0.5880, 0    ),
	c(-0.0612, 0    , -0.0612, 0    ),
	c( 0.3239, 0    ,  0.3239, 0    ),
	c( 0.6936, 0    ,  0.6936, 0    ),
	c( 0.8856, 0    ,  0.8856, 0    ),
	c( 1.0995, 0    ,  1.0995, 0    ),
	c( 1.3637, 0    ,  1.3637, 0    ),
	c( 1.5031, 0    ,  1.5031, 0    ),
	c( 1.7498, 0    ,  1.7498, 0    ),
	c( 2.0733, 0    ,  2.0733, 0    ),
	c( 2.3768, 0    ,  2.3768, 0    ))

Mx1R <- rbind(
	c(1.0000,  0.2955,  0.1268, 0.0760),
	c(0.2955,  1.0000, -0.0011, 0.1869),
	c(0.1268, -0.0011,  1.0000, 0.4377),
	c(0.0760,  0.1869,  0.4377, 1.0000))

nameList <- names(data)
# Define the model
model <- mxModel()
model <- mxModel(model, mxMatrix("Stand", name = "R", # values=c(.2955, .1268, -.0011, .0760, .1869, .4377), 
																 nrow = nvar, ncol = nvar, free=TRUE))
model <- mxModel(model, mxMatrix("Zero", name = "M", nrow = 1, ncol = nvar, free=FALSE))
model <- mxModel(model, mxMatrix("Full", 
																 name="thresh", 
																 # values = Mx1Threshold,
																 values=cbind(
																 	seq(-1.9, 1.9, length.out=nthresh2),          # t1Neur1: 12 thresholds evenly spaced from -1.9 to 1.9
																 	c(rep(1, nthresh1), rep(0, diff)),               # t1mddd4l: 1 threshold at 1
																 	seq(-1.9, 1.9, length.out=nthresh2),          # t2Neur1: 12 thresholds same as t1Neur1
																 	c(rep(1, nthresh1), rep(0, diff))                # t2mddd4l: 1 threshold same as t1mddd4l
																 ),
																 free = c(rep(c( rep(TRUE, nthresh2), 
																 								rep(TRUE, nthresh1), rep(FALSE, diff)
																 ), 2)), 
																 labels = rep(c(paste("neur", 1:nthresh2, sep=""),
																 							 paste("mddd4l", 1:nthresh1, sep=""), rep(NA, diff))
																 )))

# Define the objective function
objective <- mxExpectationNormal(covariance="R", means="M", dimnames=nameList, thresholds="thresh")

# Define the observed covariance matrix
dataMatrix <- mxData(data, type='raw')

# Add the objective function and the data to the model
model <- mxModel(model, objective, dataMatrix, mxFitFunctionML())

# Run the job
modelOut <- mxRun(model)
summary(modelOut)

######################### Nelder-Mead stuff:

#First make sure all four methods of simplex initialization run 
#(providing a matrix for the initial simplex is tested in another script):

plan <- omxDefaultComputePlan()
plan$steps$GD <- mxComputeNelderMead(xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",
																		 iniSimplexEdge=0.5, doPseudoHessian=T)
m2 <- mxModel(model,plan)
m2o <- mxRun(m2)
summary(m2o)
m2o$compute$steps$GD$output$finalFitValues
m2o$output$iterations

plan$steps$GD <- mxComputeNelderMead(xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="right",
																		 iniSimplexEdge=0.5, doPseudoHessian=T)
m3 <- mxModel(model,plan)
m3o <- mxRun(m3)
summary(m3o)
#^^^Not as good as m2
m3o$output$iterations

plan$steps$GD <- mxComputeNelderMead(xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,
																		 iniSimplexType="smartRight",iniSimplexEdge=0.5, doPseudoHessian=T)
m4 <- mxModel(model,plan)
m4o <- mxRun(m4)
summary(m4o)
m4o$output$iterations

set.seed(170301)
plan$steps$GD <- mxComputeNelderMead(xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="random",
																		 iniSimplexEdge=0.5, doPseudoHessian=T)
m5 <- mxModel(model,plan)
m5o <- mxRun(m5)
summary(m5o)
#^^^Not very good fit
m5o$output$iterations


#Test greedyMinimize:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	greedyMinimize=TRUE, doPseudoHessian=T)
m6 <- mxModel(model,plan)
m6o <- mxRun(m6)
summary(m6o) #<--Nice
m6o$output$iterations
#On AMD64 Linux/GNU, this solution passes the tests from models/passing/OrdinalTest.R :
omxCheckCloseEnough(mxEval(thresh, m6o)[,1], Mx1Threshold[,1], 0.03)
omxCheckCloseEnough(mxEval(thresh, m6o)[1,2], Mx1Threshold[1,2], 0.01)
omxCheckCloseEnough(mxEval(R, m6o), Mx1R, 0.01)
omxCheckCloseEnough(m6o$output$Minus2LogLikelihood, 4081.48, 0.08)
omxCheckCloseEnough(
	sqrt(diag(chol2inv(chol(m6o$compute$steps$GD$output$pseudoHessian)))),
	as.vector(m6o$output$standardErrors),
	5e-3
)
omxCheckTrue(all(eigen(m6o$output$hessian,T,T)$values>0))
omxCheckTrue(all(eigen(m6o$compute$steps$GD$output$pseudoHessian,T,T)$values>0))


#Test altContraction:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	altContraction=TRUE, doPseudoHessian=T)
m7 <- mxModel(model,plan)
m7o <- mxRun(m7)
summary(m7o)
m7o$output$iterations


#Test degenLimit:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	degenLimit=pi/180, doPseudoHessian=T)
m8 <- mxModel(model,plan)
m8o <- mxRun(m8)
summary(m8o)
m8o$output$iterations


#Test stagnCtrl:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	stagnCtrl=c(10,10), doPseudoHessian=T)
m9 <- mxModel(model,plan)
m9o <- mxRun(m9)
summary(m9o)
m9o$output$iterations


#Try turning off validation restart:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	validationRestart=FALSE, doPseudoHessian=T)
m10 <- mxModel(model,plan)
m10o <- mxRun(m10)
summary(m10o)
m10o$output$iterations


#Make sure the model runs when changing the default values of transformation coefficients:

#alpha:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, alpha=0.9, 
	doPseudoHessian=T)
m11 <- mxModel(model,plan)
m11o <- mxRun(m11)
summary(m11o)
m11o$output$iterations #<--maxed out

#betao:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	betao=0.4, doPseudoHessian=T)
m12 <- mxModel(model,plan)
m12o <- mxRun(m12)
summary(m12o)
m12o$output$iterations #<--maxed out

#betai:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	betai=0.4, doPseudoHessian=T)
m13 <- mxModel(model,plan)
m13o <- mxRun(m13)
summary(m13o)
m13o$output$iterations

#gamma:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	gamma=1.5, doPseudoHessian=T)
m14 <- mxModel(model,plan)
m14o <- mxRun(m14)
summary(m14o) #<--Nice
m14o$output$iterations
#On AMD64 Linux/GNU, this solution passes the tests from models/passing/OrdinalTest.R :
omxCheckCloseEnough(mxEval(thresh, m14o)[,1], Mx1Threshold[,1], 0.03)
omxCheckCloseEnough(mxEval(thresh, m14o)[1,2], Mx1Threshold[1,2], 0.01)
omxCheckCloseEnough(mxEval(R, m14o), Mx1R, 0.01)
omxCheckCloseEnough(m14o$output$Minus2LogLikelihood, 4081.48, 0.08)
omxCheckCloseEnough(
	sqrt(diag(chol2inv(chol(m14o$compute$steps$GD$output$pseudoHessian)))),
	as.vector(m14o$output$standardErrors),
	5e-3
)
omxCheckTrue(all(eigen(m14o$output$hessian,T,T)$values>0))
omxCheckTrue(all(eigen(m14o$compute$steps$GD$output$pseudoHessian,T,T)$values>0))

#gamma<=0:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	gamma=-1.5, doPseudoHessian=T)
m15 <- mxModel(model,plan)
m15o <- mxRun(m15)
summary(m15o) #<--Nice
m15o$output$iterations
#On AMD64 Linux/GNU, this solution passes the tests from models/passing/OrdinalTest.R :
omxCheckCloseEnough(mxEval(thresh, m15o)[,1], Mx1Threshold[,1], 0.03)
omxCheckCloseEnough(mxEval(thresh, m15o)[1,2], Mx1Threshold[1,2], 0.01)
omxCheckCloseEnough(mxEval(R, m15o), Mx1R, 0.01)
omxCheckCloseEnough(m15o$output$Minus2LogLikelihood, 4081.48, 0.08)
omxCheckCloseEnough(
	sqrt(diag(chol2inv(chol(m15o$compute$steps$GD$output$pseudoHessian)))),
	as.vector(m15o$output$standardErrors),
	5e-3
)
omxCheckTrue(all(eigen(m15o$output$hessian,T,T)$values>0))
omxCheckTrue(all(eigen(m15o$compute$steps$GD$output$pseudoHessian,T,T)$values>0))


#sigma:
plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	sigma=0.4, doPseudoHessian=T)
m16 <- mxModel(model,plan)
m16o <- mxRun(m16)
summary(m16o) 
m16o$output$iterations

#sigma<=0:
#Under 32-bit Windows, if using a non-random simplex for m17, 
#Nelder-Mead gets stuck in a loop of restarting the simplex
#to the same state over and over, every iteration 
#(although using a random simplex doesn't yield a good solution, at least not with this script's RNG seed).
#All the literature I've read says that shrink transformations are rare, but that has not been my experience
#so far; the user turns off shrinks at his/her own peril:
if(.Platform$OS.type=="windows" && .Platform$r_arch=="i386"){
	plan$steps$GD <- mxComputeNelderMead(
		xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="random",iniSimplexEdge=0.5, 
		sigma=-0.4, doPseudoHessian=T)
} else{
	plan$steps$GD <- mxComputeNelderMead(
		xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
		sigma=-0.4, doPseudoHessian=T)	
}
m17 <- mxModel(model,plan)
m17o <- mxRun(m17)
summary(m17o) 
m17o$output$iterations

plan$steps$GD <- mxComputeNelderMead(
	xTolProx=1e-12,fTolProx=1e-8,maxIter=10000L,iniSimplexType="regular",iniSimplexEdge=0.5, 
	doPseudoHessian=T, centerIniSimplex=TRUE)
m18 <- mxModel(model,plan)
m18o <- mxRun(m18)
summary(m18o) 
m18o$output$iterations
