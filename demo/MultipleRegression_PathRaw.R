require(OpenMx)
myRegDataRaw<-data("myRegData.txt")

myRegDataRaw<-myRegDataRaw[,c("x","y","z")]

multiRegModel <- mxModel("Multiple Regression -- Path Specification", 
      type="RAM",
      mxData(
          observed=myRegDataRaw, 
          type="raw"
      ),
      manifestVars=c("x", "y", "z"),
      # variance paths
      mxPath(
          from=c("x", "y", "z"), 
          arrows=2,
          free=TRUE, 
          values = c(1, 1, 1),
          labels=c("varx", "residual", "varz")
      ),
      # covariance of x and z
      mxPath(
          from="x",
          to="z",
          arrows=2,
          free=TRUE,
          values=0.5,
          labels="covxz"
      ), 
      # regression weights
      mxPath(
          from=c("x","z"),
          to="y",
          arrows=1,
          free=TRUE,
          values=1,
          labels=c("betax","betaz")
      ), 
      # means and intercepts
      mxPath(
          from="one",
          to=c("x", "y", "z"),
          arrows=1,
          free=TRUE,
          values=c(1, 1),
          labels=c("meanx", "beta0", "meanz")
      )
  ) # close model
      
multiRegFit<-mxRun(multiRegModel)

multiRegFit@output

# Old Mx Output
omxCheckCloseEnough(multiRegFit@output$estimate[["beta0"]], 1.6332, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["betax"]], 0.4246, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["betaz"]], 0.2260, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["residual"]], 0.6267, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["varx"]], 1.1053, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["varz"]], 0.8275, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["covxz"]], 0.2862, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(multiRegFit@output$estimate[["meanz"]], 4.0611, 0.001)

# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["beta0"]], 1.6331, 0.001)
# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["betax"]], 0.4246, 0.001)
# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["betaz"]], 0.2260, 0.001)
# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["residual"]], 0.646, 0.001)
# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["varx"]], 1.116, 0.001) 
# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["varz"]], 0.836, 0.001)
# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["covxz"]], 0.289, 0.001)
# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["meanx"]], 0.054, 0.001)
# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["meanz"]], 4.061, 0.001)
