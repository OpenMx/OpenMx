require(OpenMx)

myRegDataRaw <- read.table("myRegData.txt",header=TRUE)

multivariateRegModel <- mxModel("MultiVariate Regression -- Path Specification", 
    type="RAM",
    mxData(
        data=myRegDataRaw, 
        type="raw"
    ),
    manifestVars=c("w", "x", "y", "z"),
    # variance paths
    mxPath(
        from=c("w", "x", "y", "z"), 
        arrows=2,
        free=TRUE, 
        values = c(1, 1, 1),
        labels=c("residualw", "varx", "residualy", "varz")
    ),
    # covariance of x and z
    mxPath(
        from="x",
        to="y",
        arrows=2,
        free=TRUE,
        values=0.5,
        labels="covxz"
    ), 
    # regression weights for y
    mxPath(
        from=c("x","z"),
        to="y",
        arrows=1,
        free=TRUE,
        values=1,
        labels=c("betayx","betayz")
    ), 
    # regression weights for w
    mxPath(
        from=c("x","z"),
        to="w",
        arrows=1,
        free=TRUE,
        values=1,
        labels=c("betawx","betawz")
    ), 
    # means and intercepts
    mxPath(
        from="one",
        to=c("w", "x", "y", "z"),
        arrows=1,
        free=TRUE,
        values=c(1, 1),
        labels=c("betaw", "meanx", "betay", "meanz")
    )
) # close model
  
multivariateRegFit <- mxRun(multivariateRegModel)

multivariateRegFit@output
  
mxSummary(multivariateRegFit)  

omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betay"]], 1.6331, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betayx"]], 0.4246, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betayz"]], 0.2260, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["residualy"]], 0.646, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betaw"]], 0.51391, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betawx"]], -0.23102, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["betawz"]], 0.51223, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["residualw"]], 0.60964, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["varx"]], 1.116, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["varz"]], 0.836, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["covxz"]], 0.289, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["meanx"]], 0.054, 0.001)
omxCheckCloseEnough(multiRegPathRaw@output$estimate[["meanz"]], 4.061, 0.001)