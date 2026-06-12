library(OpenMx)
library(testthat)

jointModel <- mxModel("parent", 
	mxModel("child", type="RAM", manifestVars="X", 
	      mxData(data.frame(X=c(1,2,1,2)), "raw"),
	      mxPath("X", arrows=2, values=1, label="B", free=TRUE),
	      mxPath("one", to="X", values=1, free=TRUE)
	),
	mxFitFunctionMultigroup("child"),
	mxMatrix("Full",1,1, values = 1, free = TRUE, name = "l", labels="lambda"),
	mxPenaltyLASSO("B", "lasso", hyperparams = "lambda")
)
fit <- expect_warning(mxPenaltySearch(jointModel))
detail <- fit$compute$steps$PS$output$detail
expect_true(is.data.frame(detail))
expect_equal(nrow(detail), 41)
expect_true("EBIC" %in% colnames(detail))
expect_true("B" %in% colnames(detail))
