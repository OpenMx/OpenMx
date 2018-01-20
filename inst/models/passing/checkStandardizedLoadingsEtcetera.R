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
data(demoOneFactor)

Rff <- mxMatrix(type="Stand",nrow=1,ncol=1,free=F,name="Rff")
L <- mxMatrix(type="Full",nrow=5,ncol=1,free=T,values=0.2,labels=paste("l",1:5,sep=""),lbound=-2,ubound=2,name="L")
I <- mxMatrix(type="Iden",nrow=5,ncol=5,name="I")
C <- mxAlgebra(L %*% Rff %*% t(L),name="C")
U <- mxAlgebra(I-(I*C),name="U")
SD <- mxMatrix(type="Full",nrow=5,ncol=1,free=T,values=0.6,
               labels=c("sd1","sd2","sd3","sd4","sd5"),lbound=0,ubound=2,name="SD")
Sigma <- mxAlgebra( (C+U) * (SD%*%t(SD)), name="Sigma", dimnames=list(colnames(demoOneFactor),colnames(demoOneFactor)))

factorModelMatrix <- mxModel(
  "OneFactorMatrix", Rff, L, I, C, U, SD, Sigma, 
  mxExpectationNormal(covariance="Sigma"),
  mxFitFunctionML(),
  mxCI("L"),
  mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
)
factorMatrixFit <- mxRun(factorModelMatrix, suppressWarnings = T)
omxCheckError(mxStandardizeRAMpaths(factorMatrixFit),"model 'OneFactorMatrix' does not use RAM expectation")

manifests <- names(demoOneFactor)
latents <- c("G")
factorModelPath <- mxModel("OneFactorPath",
                           type="RAM",
                           manifestVars = manifests,
                           latentVars = latents,
                           mxPath(from=latents, to=manifests,
                                  labels=paste("L",1:5,sep="")),
                           mxPath(from=manifests, arrows=2),
                           mxPath(from=latents, arrows=2,
                                  free=FALSE, values=1.0),
                           mxData(cov(demoOneFactor), type="cov",
                                  numObs=500))
omxCheckWarning(
	mxStandardizeRAMpaths(factorModelPath,T),
  "standard errors will not be computed because model 'OneFactorPath' has not yet been run, and no matrix was provided for argument 'cov'")
factorPathFit <- mxRun(factorModelPath, suppressWarnings = T)

( matrixFitPar <- summary(factorMatrixFit)$parameters )
( pathFitPar <- summary(factorPathFit)$parameters )
( zpath <- mxStandardizeRAMpaths(factorPathFit,T) )

if (0) {
  options(digits = 12)
  print(factorPathFit$output$fit)
}
omxCheckCloseEnough(factorPathFit$output$fit, -3660.59669, 1e-3)
if (0) {
  max(abs(matrixFitPar$Estimate[1:5] - zpath$Std.Value[1:5]))
}
omxCheckCloseEnough(matrixFitPar$Estimate[1:5],zpath$Std.Value[1:5],5e-5)
omxCheckCloseEnough(matrixFitPar$Std.Error[1:5],zpath$Std.SE[1:5],5e-5)
omxCheckCloseEnough(
  zpath$Std.Value,
  c(0.891309322741553,0.932554560675801,0.943846635668634,0.962362487593728,0.972555617431608,0.205567691193994,
    0.130341991362763,0.109153528337001,0.0738584424724116,0.0541355710022234,1),
  3e-4)
omxCheckCloseEnough(
  zpath$Std.SE,
  c(0.0097233354792961,0.00640855820968981,0.00551439004191031,0.00400750800616116,0.00327547260416904,
    0.0173329991211512,0.0119526603725942,0.0104094769784555,0.00771335074796404,0.00637115856365944,0),
  5e-5)


#Check errors and warnings:
pointlessAlg <- mxModel(factorPathFit, mxAlgebra(1))
omxCheckWarning(
	mxStandardizeRAMpaths(pointlessAlg),
	"MxModel 'OneFactorPath' was modified since it was run."
)
pointlessConstraint <- mxModel(factorModelPath, mxConstraint(1==1))
pointlessConstraint <- mxRun(pointlessConstraint)
omxCheckWarning(
	mxStandardizeRAMpaths(pointlessConstraint,T),
	"standard errors will not be computed because model 'OneFactorPath' contains at least one mxConstraint"
)
#NPSOL populates the 'hessian' slot of the output with its own final Hessian when there is no MxComputeNumericDeriv step:
if(mxOption(NULL,"Default optimizer") != "NPSOL"){
	plan <- omxDefaultComputePlan()
	plan$steps <- list(plan$steps$GD, plan$steps$RE)
	nohess <- mxModel(factorModelPath, plan)
	nohess <- mxRun(nohess)
	omxCheckWarning(
		mxStandardizeRAMpaths(nohess,T),
		"argument 'SE=TRUE' requires model to have a nonempty 'hessian' output slot, or a non-NULL value for argument 'cov'; continuing with 'SE' coerced to 'FALSE'"
	)
}

#Make more models and check mxStandardizeRAMpaths()'s output for multigroup:
data("twinData", package="OpenMx")
selVars <- c('bmi1','bmi2')
aceVars <- c("A1","C1","E1","A2","C2","E2")
mzfData <- as.matrix(subset(twinData, zyg==1, selVars))
dzfData <- as.matrix(subset(twinData, zyg==3, selVars))
cov(mzfData, use="pairwise.complete.obs")
cov(dzfData, use="pairwise.complete.obs")


share <- mxModel("share", type="RAM",
                 manifestVars=selVars,
                 latentVars=aceVars,
                 mxPath(from=aceVars, arrows=2, free=FALSE, values=1),
                 mxPath(from="one", to=aceVars, arrows=1, free=FALSE, values=0),
                 mxPath(from="one", to=selVars, arrows=1, free=TRUE, values=20, labels= c("mean","mean")),
                 mxPath(from="C1",  to="C2",    arrows=2, free=FALSE, values=1),
                 mxPath(from=c("A1","C1","E1"), to="bmi1", arrows=1, free=TRUE, values=.6, label=c("a","c","e")),
                 mxPath(from=c("A2","C2","E2"), to="bmi2", arrows=1, free=TRUE, values=.6, label=c("a","c","e"))
)

MZ = mxModel(share, name="MZ",
             mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=1),
             mxData(mzfData, type="raw")
)

DZ = mxModel(share, name="DZ",
             mxPath(from="A1", to="A2", arrows=2, free=FALSE, values=.5),
             mxData(dzfData, type="raw") 
)

twinACEModel <- mxModel("twinACE", MZ, DZ,
                        mxAlgebra(MZ.objective + DZ.objective, name="twin"), 
                        mxFitFunctionAlgebra("twin")
)

myLongitudinalDataCov <- matrix(
  c(6.362, 4.344, 4.915,  5.045,  5.966,
    4.344, 7.241, 5.825,  6.181,  7.252,
    4.915, 5.825, 9.348,  7.727,  8.968,
    5.045, 6.181, 7.727, 10.821, 10.135,
    5.966, 7.252, 8.968, 10.135, 14.220),
  nrow=5,
  dimnames=list(
    c("x1","x2","x3","x4","x5"),
    c("x1","x2","x3","x4","x5"))
)

myLongitudinalDataMean <- c(9.864, 11.812, 13.612, 15.317, 17.178)
names(myLongitudinalDataMean) <- c("x1","x2","x3","x4","x5")

growthCurveModel <- mxModel("LinearGrowthCurveModel_MatrixSpecification", 
                            mxData(myLongitudinalDataCov, 
                                   type="cov", 
                                   numObs=500,
                                   mean=myLongitudinalDataMean),
                            mxMatrix(
                              type="Full",
                              nrow=7, 
                              ncol=7,
                              free=F,
                              values=c(0,0,0,0,0,1,0,
                                       0,0,0,0,0,1,1,
                                       0,0,0,0,0,1,2,
                                       0,0,0,0,0,1,3,
                                       0,0,0,0,0,1,4,
                                       0,0,0,0,0,0,0,
                                       0,0,0,0,0,0,0),
                              byrow=TRUE,
                              name="AA"),
                            mxMatrix(
                              type="Symm",
                              nrow=7,
                              ncol=7,
                              free=c(T, F, F, F, F, F, F,
                                     F, T, F, F, F, F, F,
                                     F, F, T, F, F, F, F,
                                     F, F, F, T, F, F, F,
                                     F, F, F, F, T, F, F,
                                     F, F, F, F, F, T, T,
                                     F, F, F, F, F, T, T),
                              values=c(0,0,0,0,0,  0,  0,
                                       0,0,0,0,0,  0,  0,
                                       0,0,0,0,0,  0,  0,
                                       0,0,0,0,0,  0,  0,
                                       0,0,0,0,0,  0,  0,
                                       0,0,0,0,0,  1,0.5,
                                       0,0,0,0,0,0.5,  1),
                              labels=c("residual", NA, NA, NA, NA, NA, NA,
                                       NA, "residual", NA, NA, NA, NA, NA,
                                       NA, NA, "residual", NA, NA, NA, NA,
                                       NA, NA, NA, "residual", NA, NA, NA,
                                       NA, NA, NA, NA, "residual", NA, NA,
                                       NA, NA, NA, NA, NA, "vari", "cov",
                                       NA, NA, NA, NA, NA, "cov", "vars"),
                              byrow= TRUE,
                              name="SS"),
                            mxMatrix(
                              type="Full",
                              nrow=5,
                              ncol=7,
                              free=F,
                              values=c(1,0,0,0,0,0,0,
                                       0,1,0,0,0,0,0,
                                       0,0,1,0,0,0,0,
                                       0,0,0,1,0,0,0,
                                       0,0,0,0,1,0,0),
                              byrow=T,
                              name="F"),
                            mxMatrix("Full",  nrow=1, ncol=7,
                                     values=c(0,0,0,0,0,1,1),
                                     free=c(F,F,F,F,F,T,T),
                                     labels=c(NA,NA,NA,NA,NA,"meani","means"),
                                     name="M"),
                            mxFitFunctionML(),mxExpectationRAM("AA","SS","F","M",dimnames=c("x1","x2","x3","x4","x5","intercept","slope"))
)

bigmod <- mxModel(
  model="foo",factorModelPath,factorModelMatrix,twinACEModel,growthCurveModel,
  mxFitFunctionMultigroup(c("OneFactorPath.fitfunction","OneFactorMatrix.fitfunction","twinACE.fitfunction",
                            "LinearGrowthCurveModel_MatrixSpecification.fitfunction"))
)
bigrun <- mxRun(bigmod, suppressWarnings = T)
zpath2 <- mxStandardizeRAMpaths(bigrun,T)
omxCheckEquals(names(zpath2),c("OneFactorPath","twinACE","LinearGrowthCurveModel_MatrixSpecification"))
omxCheckEquals(names(zpath2[[2]]),c("MZ","DZ"))


bigmod2 <- mxModel(
  model="foo",mxModel(factorModelPath,independent=T),mxModel(factorModelMatrix,independent=T),
  twinACEModel,growthCurveModel,
  mxFitFunctionMultigroup(c("twinACE.fitfunction","LinearGrowthCurveModel_MatrixSpecification.fitfunction"))
)
bigrun2 <- mxRun(bigmod2, suppressWarnings = T)
zpath3 <- mxStandardizeRAMpaths(bigrun2,T)
omxCheckEquals(names(zpath3),c("OneFactorPath","twinACE","LinearGrowthCurveModel_MatrixSpecification"))
omxCheckEquals(names(zpath3[[2]]),c("MZ","DZ"))
omxCheckEquals(sum(zpath3$OneFactorPath==zpath,na.rm=T),93)


bigmod3 <- mxModel(
  model="foo",mxModel(factorModelPath,independent=T),mxModel(factorModelMatrix,independent=T),
  mxModel(twinACEModel,independent=T),mxModel(growthCurveModel,independent=T))
bigrun3 <- mxRun(bigmod3, suppressWarnings = T)
zpath4 <- mxStandardizeRAMpaths(bigrun3,T)
omxCheckEquals(sum(zpath4$OneFactorPath==zpath3$OneFactorPath,na.rm=T),93)
omxCheckEquals(names(zpath4),c("OneFactorPath","twinACE","LinearGrowthCurveModel_MatrixSpecification"))
omxCheckEquals(names(zpath4[[2]]),c("MZ","DZ"))

#Regression test to ensure (1) that if the "container" model uses RAM expectation, its results are included in the standardized-paths output,
#and (2) that its element of the output list is named after it:
littlemod <- mxModel(factorModelPath, mxModel(growthCurveModel, independent=T))
littlerun <- mxRun(littlemod)
zpath5 <- mxStandardizeRAMpaths(littlerun)
omxCheckEquals(names(zpath5)[1],factorModelPath$name)
omxCheckEquals(zpath4$OneFactorPath$Std.Value,zpath5$OneFactorPath$Std.Value)
omxCheckEquals(zpath4$LinearGrowthCurveModel_MatrixSpecification$Std.Value,zpath5$LinearGrowthCurveModel_MatrixSpecification$Std.Value)

# #Test mxTryHard(), since the initial run of 'bigmod' gets Code Red:
# set.seed(420)
# bigrun_again <- mxTryHard(bigmod)
# omxCheckEquals(bigrun_again$output$status$code,0)
