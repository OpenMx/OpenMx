require(OpenMx)

myRegDataCov <- matrix(
    c(0.808,-0.110, 0.089, 0.361,
     -0.110, 1.116, 0.539, 0.289,
      0.089, 0.539, 0.933, 0.312,
      0.361, 0.289, 0.312, 0.836),
    nrow=4,
    dimnames=list(
      c("w","x","y","z"),
      c("w","x","y","z"))
)
 
myRegDataMeans <- c(2.582, 0.054, 2.574, 4.061)

MultipleDataCov <- myRegDataCov[c("x","y","z"),c("x","y","z")]	
MultipleDataMeans <- myRegDataMeans[c(2,3,4)]

multiRegModel <- mxModel("Multiple Regression -- Path Specification", 
      type="RAM",
      mxData(
          observed=MultipleDataCov, 
          type="cov",
          numObs=100,
          means=MultipleDataMeans
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
      
multipleRegPathCov<-mxRun(multiRegModel)

multipleRegPathCov@output

# Old Mx Output
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["beta0"]], 1.6312, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["betax"]], 0.4243, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["betaz"]], 0.2265, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["residual"]], 0.6336, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["varx"]], 1.1160, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["varz"]], 0.8360, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["covxz"]], 0.2890, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["meanx"]], 0.0540, 0.001)
omxCheckCloseEnough(multipleRegPathCov@output$estimate[["meanz"]], 4.0610, 0.001)

# omxCheckCloseEnough(multipleRegPathCov@output$estimate[["beta0"]], 1.6331, 0.001)
# omxCheckCloseEnough(multipleRegPathCov@output$estimate[["betax"]], 0.4246, 0.001)
# omxCheckCloseEnough(multipleRegPathCov@output$estimate[["betaz"]], 0.2260, 0.001)
# omxCheckCloseEnough(multipleRegPathCov@output$estimate[["residual"]], 0.646, 0.001)
# omxCheckCloseEnough(multipleRegPathCov@output$estimate[["varx"]], 1.116, 0.001)
# omxCheckCloseEnough(multipleRegPathCov@output$estimate[["varz"]], 0.836, 0.001)
# omxCheckCloseEnough(multipleRegPathCov@output$estimate[["covxz"]], 0.289, 0.001)
# omxCheckCloseEnough(multipleRegPathCov@output$estimate[["meanx"]], 0.054, 0.001)
# omxCheckCloseEnough(multipleRegPathRaw@output$estimate[["meanz"]], 4.061, 0.001)
