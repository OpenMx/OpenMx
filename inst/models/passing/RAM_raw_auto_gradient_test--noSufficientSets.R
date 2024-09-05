require(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")

mydat <- demoOneFactor[1:25,]
mydat[1,1] <- NA
mydat[2,2] <- NA
mydat[3,3] <- NA
mydat[4,4] <- NA
mydat[5,5] <- NA
mydat[6,c(1,2)] <- NA
mydat[7,c(1,3)] <- NA
mydat[8,c(1,4)] <- NA
mydat[9,c(1,5)] <- NA
mydat[10,c(2,3)] <- NA
mydat[11,c(2,4)] <- NA
mydat[12,c(2,5)] <- NA
mydat[13,c(3,4)] <- NA
mydat[14,c(3,5)] <- NA
mydat[15,c(4,5)] <- NA
mydat[16,c(1,2,3)] <- NA
mydat[17,c(1,2,4)] <- NA
mydat[18,c(1,2,5)] <- NA
mydat[19,c(1,3,4)] <- NA
mydat[20,c(1,3,5)] <- NA
mydat[21,c(1,4,5)] <- NA
mydat[22,c(2,3,4)] <- NA
mydat[23,c(2,3,5)] <- NA
mydat[24,c(2,4,5)] <- NA
mydat[25,c(3,4,5)] <- NA
plan <- mxComputeSequence(list(
	GD=mxComputeGradientDescent(verbose=5L),ND=mxComputeNumericDeriv(verbose=5L),SE=mxComputeStandardError(),
	HQ=mxComputeHessianQuality(verbose=5L),RD=mxComputeReportDeriv(),RE=mxComputeReportExpectation()
))
factorModel <- mxModel(
	"One Factor",
	type="RAM",
	#plan,
	manifestVars = manifests,
	latentVars = latents,
	mxPath(from=latents, to=manifests,values=0.8,labels="l"),
	mxPath(from=manifests, arrows=2,values=1,labels="u"),
	mxPath(from=latents, arrows=2,
				 free=FALSE, values=1.0),
	mxPath(from="one",to=manifests,values=0.1,labels="m"),
	mxData(mydat, type="raw")
)
mxOption(NULL,"Analytic gradients","Yes")
fmf5 <- mxRun(factorModel)
summary(fmf5)
mxOption(NULL,"Analytic gradients","No")
fmf6 <- mxRun(factorModel)
summary(fmf6)
omxCheckCloseEnough(coef(fmf5),coef(fmf6),5e-6)
omxCheckCloseEnough(fmf5$output$standardErrors,fmf6$output$standardErrors,5e-7)
omxCheckCloseEnough(fmf5$output$fit,fmf6$output$fit,1e-8)
omxCheckTrue(fmf5$output$status$code==0)
#Check for zero gradient:
omxCheckCloseEnough(fmf5$output$gradient,c(0,0,0),5e-4)
omxCheckCloseEnough(fmf6$output$gradient,c(0,0,0),5e-4)
#Using analytic derivatives should be faster:
omxCheckTrue(fmf5$output$iterations <= fmf6$output$iterations)
omxCheckTrue(fmf5$output$evaluations < fmf6$output$evaluations)
if(0){
	omxCheckTrue(summary(fmf5)$wallTime < summary(fmf6)$wallTime)
}
