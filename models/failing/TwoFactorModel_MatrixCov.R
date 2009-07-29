require(OpenMx)

myFADataCov <- matrix(
      c(0.997, 0.642, 0.611, 0.672, 0.637, 0.677, 0.342, 0.299, 0.337,
        0.642, 1.025, 0.608, 0.668, 0.643, 0.676, 0.273, 0.282, 0.287,
        0.611, 0.608, 0.984, 0.633, 0.657, 0.626, 0.286, 0.287, 0.264,
        0.672, 0.668, 0.633, 1.003, 0.676, 0.665, 0.330, 0.290, 0.274,
        0.637, 0.643, 0.657, 0.676, 1.028, 0.654, 0.328, 0.317, 0.331,
        0.677, 0.676, 0.626, 0.665, 0.654, 1.020, 0.323, 0.341, 0.349,
        0.342, 0.273, 0.286, 0.330, 0.328, 0.323, 0.993, 0.472, 0.467,
        0.299, 0.282, 0.287, 0.290, 0.317, 0.341, 0.472, 0.978, 0.507,
        0.337, 0.287, 0.264, 0.274, 0.331, 0.349, 0.467, 0.507, 1.059),
      nrow=9,
      dimnames=list(
          c("x1", "x2", "x3", "x4", "x5", "x6", "y1", "y2", "y3"),
          c("x1", "x2", "x3", "x4", "x5", "x6", "y1", "y2", "y3")),
      )

twoFactorCov <- myFADataCov[c("x1", "x2", "x3", "y1", "y2", "y3"),c("x1", "x2", "x3", "y1", "y2", "y3")]
  
myFADataMeans <- c(2.988, 3.011, 2.986, 3.053, 3.016, 3.010, 2.955, 2.956, 2.967)
  
twoFactorMeans <- myFADataMeans[c(1:3,7:9)]

twoFactorModel <- mxModel("Two Factor Model - Matrix Specification", 
    mxData(
        data=twoFactorCov, 
        type="cov", 
        numObs=500,
        means=twoFactorMeans
        ),
    mxMatrix("Full", nrow=8, ncol=8,
        values=c(0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,1,
                 0,0,0,0,0,0,0,0,
                 0,0,0,0,0,0,0,0),
        free=c(F, F, F, F, F, F, F, F,
               F, F, F, F, F, F, T, F,
               F, F, F, F, F, F, T, F,
               F, F, F, F, F, F, F, F,
               F, F, F, F, F, F, F, T,
               F, F, F, F, F, F, F, T,
               F, F, F, F, F, F, F, F,
               F, F, F, F, F, F, F, F),
        labels=c(NA,NA,NA,NA,NA,NA,"l1", NA,
                 NA,NA,NA,NA,NA,NA,"l2", NA,
                 NA,NA,NA,NA,NA,NA,"l3", NA,
                 NA,NA,NA,NA,NA,NA, NA,"l4",
                 NA,NA,NA,NA,NA,NA, NA,"l5",
                 NA,NA,NA,NA,NA,NA, NA,"l6",
                 NA,NA,NA,NA,NA,NA, NA, NA,
                 NA,NA,NA,NA,NA,NA, NA, NA),
        byrow=TRUE,
        name="A"),
    mxMatrix("Symm", nrow=8, ncol=8, 
        values=c(1,0,0,0,0,0, 0, 0,
                 0,1,0,0,0,0, 0, 0,
                 0,0,1,0,0,0, 0, 0,
                 0,0,0,1,0,0, 0, 0,
                 0,0,0,0,1,0, 0, 0,
                 0,0,0,0,0,1, 0, 0,
                 0,0,0,0,0,0, 1,.5,
                 0,0,0,0,0,0,.5, 1),
        free=c(T, F, F, F, F, F, F, F,
               F, T, F, F, F, F, F, F,
               F, F, T, F, F, F, F, F,
               F, F, F, T, F, F, F, F,
               F, F, F, F, T, F, F, F,
               F, F, F, F, F, T, F, F,
               F, F, F, F, F, F, T, T,
               F, F, F, F, F, F, T, T),
        labels=c("e1", NA,   NA,   NA,   NA,   NA,    NA,    NA,
                 NA, "e2",   NA,   NA,   NA,   NA,    NA,    NA,
                 NA,   NA, "e3",   NA,   NA,   NA,    NA,    NA,
                 NA,   NA,   NA, "e4",   NA,   NA,    NA,    NA,
                 NA,   NA,   NA,   NA, "e5",   NA,    NA,    NA,
                 NA,   NA,   NA,   NA,   NA, "e6",    NA,    NA,
                 NA,   NA,   NA,   NA,   NA,   NA, "varF1", "cov",
                 NA,   NA,   NA,   NA,   NA,   NA, "cov", "varF2"),
        byrow=TRUE,
        name="S"),
    mxMatrix("Full",  nrow=6, ncol=8,
        free=F,
        values=c(1,0,0,0,0,0,0,0,
                 0,1,0,0,0,0,0,0,
                 0,0,1,0,0,0,0,0,
                 0,0,0,1,0,0,0,0,
                 0,0,0,0,1,0,0,0,
                 0,0,0,0,0,1,0,0),
        byrow=T,
        name="F"),
    mxMatrix("Full", nrow=1, ncol=8,
        values=c(1,1,1,1,1,1,0,0),
        free=c(T,T,T,T,T,T,F,F),
        labels=c("meanx1","meanx2","meanx3",
                 "meanx4","meanx5","meanx6",
                  NA,NA),
        name="M"),
    mxRAMObjective("A","S","F","M")
)
      
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

