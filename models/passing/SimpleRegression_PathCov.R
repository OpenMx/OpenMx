myRegDataCov <- matrix(
	c(1.116, 0.539,
	  0.539, 0.933),
	nrow=2,
	dimnames=list(c("x", "y"), c("x", "y")))
	
myRegDataMeans<-c(0.05416, 2.57393)

# Create an MxModel object
uniRegModel <- mxModel("Simple Regression -- Path Specification", 
    type="RAM",
    mxData(
        data=myRegDataCov, 
        type="cov", 
        numObs=100,
        means=myRegDataMeans 
    ),
    manifestVars=c("x", "y"),
    # variances paths
    mxPath(
        from=c("x", "y"), 
        arrows=2,
        free=TRUE, 
        values = c(1, 1),
        labels=c("varx", "residual")
    ),
    mxPath(
        from="x",
        to="y",
        arrows=1,
        free=TRUE,
        values=1,
        labels="beta1"
    ), # regression weight
    mxPath(
        from="one",
        to=c("x", "y"),
        arrows=1,
        free=TRUE,
        values=c(1, 1),
        labels=c("meanx", "beta0")
    ) # means
) # close model

      
uniRegFit <- mxRun(uniRegModel)

# commented out until the summary function is implemented
# summary(uniRegFit)

uniRegFit@output

omxCheckCloseEnough(uniRegFit@output$estimate[["beta0"]], 2.54776, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["beta1"]], 0.48312, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["residual"]], 0.672, 0.01)
omxCheckCloseEnough(uniRegFit@output$estimate[["meanx"]], 0.05412, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["varx"]], 1.11654, 0.001)