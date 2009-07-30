require(OpenMx)

myFADataRaw <- read.table("myFAData.txt", header=T)

twoFactorRaw <- myFADataRaw[,c("x1","x2","x3","y1","y2","y3")]

twoFactorModel <- mxModel("Two Factor Model - Path", 
    type="RAM",
    mxData(
        observed=twoFactorRaw, 
        type="raw",
        ),
    manifestVars=c("x1", "x2", "x3", "y1", "y2", "y3"),
    latentVars=c("F1","F2"),
    # residual variances
    mxPath(from=c("x1", "x2", "x3", "y1", "y2", "y3"),
        arrows=2,
        free=TRUE,
        values=c(1,1,1,1,1,1),
        labels=c("e1","e2","e3","e4","e5","e6")
    ),
    # latent variances and covaraince
    mxPath(from=c("F1","F2"),
        arrows=2,
        all=2,
        free=TRUE,
        values=c(1, .5,
                .5, 1),
        labels=c("varF1", "cov", "cov", "varF2")
    ), 
    # factor loadings for x variables
    mxPath(from="F1",
        to=c("x1","x2","x3"),
        arrows=1,
        free=c(FALSE,TRUE,TRUE),
        values=c(1,1,1),
        labels=c("l1","l2","l3")
    ),
    # factor loadings for y variables
    mxPath(from="F2",
        to=c("y1","y2","y3"),
        arrows=1,
        free=c(FALSE,TRUE,TRUE),
        values=c(1,1,1),
        labels=c("l4","l5","l6")
    ),
    # means
    mxPath(from="one",
        to=c("x1","x2","x3","x4","x5","x6","F1","F2"),
        arrows=1,
        free=c(T ,T, T, T, T, T, F, F),
        values=c(1,1,1,1,1,1,0,0),
        labels=c("meanx1","meanx2","meanx3",
                 "meany1","meany2","meany3",
                  NA,NA)
    ) # means
) # close model
      
twoFactorFit <- mxRun(twoFactorModel)

omxCheckCloseEnough(twoFactorFit@output$estimate[["l2"]], 0.999, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l3"]], 0.959, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l4"]], 1.028, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l5"]], 1.008, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l6"]], 1.021, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["varF1"]], 0.645, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e1"]], 0.350, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e2"]], 0.379, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e3"]], 0.389, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e4"]], 0.320, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e5"]], 0.370, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e6"]], 0.346, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx1"]], 2.988, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx2"]], 3.011, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx3"]], 2.986, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx4"]], 3.053, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx5"]], 3.016, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx6"]], 3.010, 0.01)
