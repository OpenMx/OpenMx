require(OpenMx)

plan <- omxDefaultComputePlan()
omxCheckEquals(names(plan$steps),c("GD","ND","SE","HQ","RD","RE"))

plan <- omxDefaultComputePlan(useOptimizer=FALSE)
omxCheckEquals(names(plan$steps),c("CO","RE"))

plan <- omxDefaultComputePlan(intervals=TRUE)
omxCheckEquals(names(plan$steps),c("GD","CI","ND","SE","HQ","RD","RE"))

plan <- omxDefaultComputePlan(modelName="foo",intervals=TRUE)
omxCheckEquals(plan$steps$GD$fitfunction,"foo.fitfunction")
omxCheckEquals(
	a=ifelse(test=mxOption(NULL,"Default optimizer")=="SLSQP",
					 yes=plan$steps$CI$plan$plan$fitfunction,
					 no=plan$steps$CI$plan$fitfunction),
	b="foo.fitfunction")
omxCheckEquals(plan$steps$ND$fitfunction,"foo.fitfunction")

ol <- options()$mxOption
ol$"Gradient algorithm" <- "forward"
ol$"Gradient iterations" <- 2L
ol$"Gradient step size" <- 1e-3
ol$"Calculate Hessian" <- "No"
ol$"Standard Errors" <- "No"
plan <- omxDefaultComputePlan(intervals=TRUE,optionList=ol)
omxCheckEquals(names(plan$steps),c("GD","CI","RD","RE"))
omxCheckEquals(plan$steps$GD$gradientAlgo,"forward")
omxCheckEquals(plan$steps$GD$gradientIterations,2L)
omxCheckEquals(plan$steps$GD$gradientStepSize,0.001)
if(mxOption(NULL,"Default optimizer")=="SLSQP"){
	omxCheckEquals(plan$steps$CI$plan$plan$gradientAlgo,"forward")
	omxCheckEquals(plan$steps$CI$plan$plan$gradientIterations,2L)
	omxCheckEquals(plan$steps$CI$plan$plan$gradientStepSize,0.001)
} else{
	omxCheckEquals(plan$steps$CI$plan$gradientAlgo,"forward")
	omxCheckEquals(plan$steps$CI$plan$gradientIterations,2L)
	omxCheckEquals(plan$steps$CI$plan$gradientStepSize,0.001)
}
