require(OpenMx)
myRegDataCov <- matrix(
    c(1.116, 0.539,
      0.539, 0.933),
    nrow=2,
    dimnames=list(c("x","y"),c("x","y"))
)
	
myRegDataMeans <- c(0.05416, 2.57393)

uniRegModel <- mxModel("Simple Regression - Matrix Specification", 
    mxData(
      data=myRegDataCov, 
      type="cov", 
      numObs=100,
      means=myRegDataMeans
    ),
    mxMatrix(
        type="Full", 
        nrow=2, 
        ncol=2,
        free=c(FALSE, FALSE,
               TRUE,  FALSE),
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
        free=c(TRUE,  FALSE,
               FALSE, TRUE),
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
        free=c(TRUE, TRUE),
        values=c(0, 0),
        labels=c("meanx", "beta0"),
        name="M"),
    mxRAMObjective("A", "S", "F", "M")
)
      
      
uniRegFit <- mxRun(uniRegModel)

uniRegFit@output


omxCheckCloseEnough(uniRegFit@output$estimate[["beta0"]], 2.54776, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["beta1"]], 0.48312, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["residual"]], 0.672, 0.01)
omxCheckCloseEnough(uniRegFit@output$estimate[["meanx"]], 0.05412, 0.001)
omxCheckCloseEnough(uniRegFit@output$estimate[["varx"]], 1.11654, 0.001)
