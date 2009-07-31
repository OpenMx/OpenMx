require(OpenMx)
myRegDataRaw<-read.table("myRegData.txt",header=TRUE)

names(myRegDataRaw)

myRegDataRaw<-myRegDataRaw[,c("x","y")]

uniRegModel<-mxModel("Simple Regression, Path Specification", 
      type="RAM",
      mxData(myRegDataRaw, type="raw"),
      manifestVars=c("x","y"),
      mxPath(from=c("x","y"), 
            arrows=2,
            free=TRUE, 
            values = c(1,1),
            labels=c("varx","residual")), # variances
      mxPath(from="x",
            to="y",
            arrows=1,
            free=TRUE,
            values=1,
            label="beta1"), # regression weight
      mxPath(from="one",
            to=c("x","y"),
            arrows=1,
            free=TRUE,
            values=c(1,1),
            labels=c("meanx","beta0")) # means
      ) # close model
      
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
# omxCheckCloseEnough(uniRegOutput@output$estimate[["residual"]], 0.672, 0.01)
# omxCheckCloseEnough(uniRegOutput@output$estimate[["meanx"]], 0.05412, 0.001)
# omxCheckCloseEnough(uniRegOutput@output$estimate[["varx"]], 1.11654, 0.001)
