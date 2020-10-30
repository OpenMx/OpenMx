require(OpenMx)
library(testthat)

plan <- omxDefaultComputePlan()
omxCheckEquals(names(plan$steps),c("GD","ND","SE","HQ","RD","RE"))

plan <- omxDefaultComputePlan(useOptimizer=FALSE)
omxCheckEquals(names(plan$steps),c("CO","RE"))

plan <- omxDefaultComputePlan(intervals=TRUE)
omxCheckEquals(names(plan$steps),c("GD","CI","ND","SE","HQ","RD","RE"))

plan <- omxDefaultComputePlan(modelName="foo",intervals=TRUE)
omxCheckEquals(plan$steps$GD$fitfunction,"foo.fitfunction")
expect_equivalent(
  ifelse(test=mxOption(NULL,"Default optimizer")=="SLSQP",
         yes=plan$steps$CI$plan$plan$fitfunction,
         no=plan$steps$CI$plan$fitfunction),
	"foo.fitfunction")
omxCheckEquals(plan$steps$ND$fitfunction,"foo.fitfunction")
