require(OpenMx)

data(myRegDataRaw)

SimpleDataRaw<-myRegDataRaw[,c("x","y")]

uniRegModel<-mxModel("Simple Regression - Matrix Specification", 
    mxData(
        observed=SimpleDataRaw,
        type="raw"),
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
        dimnames=list(c("x","y"),c("x","y")),
        name="F"
    ),
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=1,
        free=c(T, T),
        values=c(0, 0),
        labels=c("meanx", "beta0"),
        dimnames=list(c("x","y"), NULL),
        name="M"),
    mxRAMObjective("A", "S", "F", "M")
)
      
uniRegOutput<-mxRun(uniRegModel)

uniRegOutput@output

# Old Mx Output
omxCheckCloseEnough(uniRegOutput@output$estimate[["beta0"]], 2.5478, 0.001)
omxCheckCloseEnough(uniRegOutput@output$estimate[["beta1"]], 0.4831, 0.001)
omxCheckCloseEnough(uniRegOutput@output$estimate[["residual"]], 0.6652, 0.001)
omxCheckCloseEnough(uniRegOutput@output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(uniRegOutput@output$estimate[["varx"]], 1.1053, 0.001)

# omxCheckCloseEnough(uniRegOutput@output$estimate[["beta0"]], 2.54776, 0.001)
# omxCheckCloseEnough(uniRegOutput@output$estimate[["beta1"]], 0.48312, 0.001)
# omxCheckCloseEnough(uniRegOutput@output$estimate[["residual"]], 0.672, 0.001)
# omxCheckCloseEnough(uniRegOutput@output$estimate[["meanx"]], 0.05412, 0.001)
# omxCheckCloseEnough(uniRegOutput@output$estimate[["varx"]], 1.11654, 0.001)
