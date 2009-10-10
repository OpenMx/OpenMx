# -----------------------------------------------------------------------
# Program: SimpleRegression_MatrixRaw.R  
#  Author: Ryne Estabrook
#    Date: 08 01 2009 
#
# Simple Regression model to estimate effect of independent on dependent variables
# Matrix style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
data(myRegDataRaw)

SimpleDataRaw <- myRegDataRaw[,c("x","y")]

#Create an MxModel object
# -----------------------------------------------------------------------
uniRegModel <- mxModel("Simple Regression -- Matrix Specification", 
    mxData(
        observed=SimpleDataRaw,
        type="raw"
    ),
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=2,
        free=c(F, F,
               T, F),
        values=c(0, 0,
                 1, 0),
        labels=c(NA,     NA,
                "beta1", NA),
        byrow=TRUE,
        name="A"
    ),
    mxMatrix(
        type="Symm", 
        nrow=2, 
        ncol=2, 
        values=c(1, 0,
                 0, 1),
        free=c(T, F,
               F, T),
        labels=c("varx", NA,
                  NA,    "residual"),
        byrow=TRUE,
        name="S"
    ),
    mxMatrix(
        type="Iden",  
        nrow=2, 
        ncol=2,
        name="F"
    ),
    mxMatrix(
        type="Full", 
        nrow=1, 
        ncol=2,
        free=c(T, T),
        values=c(0, 0),
        labels=c("meanx", "beta0"),
        name="M"),
    mxRAMObjective("A", "S", "F", "M")
)
      
uniRegFit<-mxRun(uniRegModel)

summary(uniRegFit)
uniRegFit@output

#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
omxCheckCloseEnough(uniRegFit@output$estimate[["beta0"]], 2.5478, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["beta1"]], 0.4831, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["residual"]], 0.6652, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["varx"]], 1.1053, 0.001)
