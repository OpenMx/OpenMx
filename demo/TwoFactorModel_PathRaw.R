# -----------------------------------------------------------------------
# Program: TwoFactorModel_PathRaw.R  
#  Author: Ryne Estabrook
#    Date: 08 01 2009 
#
# Two Factor model to estimate factor loadings, residual variances and means
# Path style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
data(myFADataRaw)

twoFactorRaw <- myFADataRaw[,c("x1","x2","x3","y1","y2","y3")]

#Create an MxModel object
# -----------------------------------------------------------------------
twoFactorModel <- mxModel("Two Factor Model -- Path Specification", 
	type="RAM",
    mxData(
    	observed=twoFactorRaw, 
    	type="raw"
    ),
    manifestVars=c("x1", "x2", "x3", "y1", "y2", "y3"),
    latentVars=c("F1","F2"),
    # residual variances
    mxPath(
    	from=c("x1", "x2", "x3", "y1", "y2", "y3"),
        arrows=2,
        free=TRUE,
        values=c(1,1,1,1,1,1),
        labels=c("e1","e2","e3","e4","e5","e6")
    ),
    # latent variances and covaraince
    mxPath(
    	from=c("F1","F2"),
        arrows=2,
        all=TRUE,
        free=TRUE,
        values=c(1, .5,
                .5, 1),
        labels=c("varF1", "cov", "cov", "varF2")
    ), 
    # factor loadings for x variables
    mxPath(
    	from="F1",
        to=c("x1","x2","x3"),
        arrows=1,
        free=c(FALSE,TRUE,TRUE),
        values=c(1,1,1),
        labels=c("l1","l2","l3")
    ),
    # factor loadings for y variables
    mxPath(
    	from="F2",
        to=c("y1","y2","y3"),
        arrows=1,
        free=c(FALSE,TRUE,TRUE),
        values=c(1,1,1),
        labels=c("l4","l5","l6")
    ),
    # means
    mxPath(
    	from="one",
        to=c("x1","x2","x3","y1","y2","y3","F1","F2"),
        arrows=1,
        free=c(T ,T, T, T, T, T, F, F),
        values=c(1,1,1,1,1,1,0,0),
        labels=c("meanx1","meanx2","meanx3",
                 "meany1","meany2","meany3",
                  NA,NA)
    )
) # close model
      
twoFactorFit <- mxRun(twoFactorModel)

summary(twoFactorFit)
twoFactorFit@output$estimate

#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
omxCheckCloseEnough(twoFactorFit@output$estimate[["l2"]], 0.9723, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l3"]], 0.9313, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l5"]], 1.0498, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l6"]], 1.0531, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["varF1"]], 0.6604, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["varF2"]], 0.4505, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["cov"]], 0.2952, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e1"]], 0.3349, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e2"]], 0.3985, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e3"]], 0.4091, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e4"]], 0.5404, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e5"]], 0.4809, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e6"]], 0.5571, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx1"]], 2.988, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx2"]], 3.0113, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx3"]], 2.9861, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meany1"]], 2.9554, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meany2"]], 2.9562, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meany3"]], 2.9673, 0.01)
