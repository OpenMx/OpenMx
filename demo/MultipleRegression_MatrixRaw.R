# -----------------------------------------------------------------------
# Program: MultipleRegression_MatrixRaw.R  
#  Author: Ryne Estabrook
#    Date: 08 01 2009 
#
# Multiple Regression model to estimate effect of independent on dependent variables
# Matrix style model input - Raw data input
#
# Revision History
#   Hermine Maes -- 10 08 2009 updated & reformatted
# -----------------------------------------------------------------------

require(OpenMx)

#Prepare Data
# -----------------------------------------------------------------------
data(myRegDataRaw)

MultipleDataRaw<-myRegDataRaw[,c("x","y","z")]

#Create an MxModel object
# -----------------------------------------------------------------------
multiRegModel <- mxModel("Multiple Regression -- Matrix Specification", 
    mxData(
    	observed=MultipleDataRaw,
    	type="raw"
   	),
    mxMatrix(
    	type="Full", 
    	nrow=3, 
    	ncol=3,
        values=c(0,0,0,
                 1,0,1,
                 0,0,0),
        free=c(F, F, F,
               T, F, T,
               F, F, F),
        labels=c(NA,     NA, NA,
                "betax", NA,"betaz",
                 NA,     NA, NA),
        byrow=TRUE,
        name = "A"
    ),
    mxMatrix(
    	type="Symm", 
    	nrow=3, 
    	ncol=3, 
        values=c(1, 0, .5,
                 0, 1, 0,
                .5, 0, 1),
        free=c(T, F, T,
               F, T, F,
               T, F, T),
        labels=c("varx",  NA,         "covxz",
                  NA,    "residual",   NA,
                 "covxz", NA,         "varz"),
        byrow=TRUE,
        name="S"
    ),
    mxMatrix(
    	type="Iden",
    	nrow=3, 
    	ncol=3,
        name="F",
    ),
    mxMatrix(
    	type="Full", 
    	nrow=1, 
    	ncol=3,
        values=c(0,0,0),
        free=c(T,T,T),
        labels=c("meanx","beta0","meanz"),
        name="M"
    ),
    mxRAMObjective("A","S","F","M")
)
      
multiRegFit<-mxRun(multiRegModel)

summary(multiRegFit)
multiRegFit@output

#Compare OpenMx results to Mx results 
# -----------------------------------------------------------------------
omxCheckCloseEnough(multiRegFit@output$estimate[["beta0"]], 1.6332, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["betax"]], 0.4246, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["betaz"]], 0.2260, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["residual"]], 0.6267, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["varx"]], 1.1053, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["varz"]], 0.8275, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["covxz"]], 0.2862, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["meanz"]], 4.0611, 0.001)
