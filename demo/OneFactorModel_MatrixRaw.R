require(OpenMx)
data(myFADataRaw)
myFADataRaw<-myFADataRaw[,c("x1","x2","x3","x4","x5","x6")]

oneFactorModel<-mxModel("Common Factor Model - Matrix Specification", 
      mxData(myFADataRaw, type="raw"),
      mxMatrix("Full", nrow=7, ncol=7,
            values=c(0,0,0,0,0,0,1,
                     0,0,0,0,0,0,1,
                     0,0,0,0,0,0,1,
                     0,0,0,0,0,0,1,
                     0,0,0,0,0,0,1,
                     0,0,0,0,0,0,1,
                     0,0,0,0,0,0,0),
            free=c(F, F, F, F, F, F, F,
                   F, F, F, F, F, F, T,
                   F, F, F, F, F, F, T,
                   F, F, F, F, F, F, T,
                   F, F, F, F, F, F, T,
                   F, F, F, F, F, F, T,
                   F, F, F, F, F, F, F),
            labels=c(NA,NA,NA,NA,NA,NA,"l1",
                     NA,NA,NA,NA,NA,NA,"l2",
                     NA,NA,NA,NA,NA,NA,"l3",
                     NA,NA,NA,NA,NA,NA,"l4",
                     NA,NA,NA,NA,NA,NA,"l5",
                     NA,NA,NA,NA,NA,NA,"l6",
                     NA,NA,NA,NA,NA,NA,NA),
            byrow=TRUE,
            name="A"),
      mxMatrix("Symm", nrow=7, ncol=7, 
            values=c(1,0,0,0,0,0,0,
                     0,1,0,0,0,0,0,
                     0,0,1,0,0,0,0,
                     0,0,0,1,0,0,0,
                     0,0,0,0,1,0,0,
                     0,0,0,0,0,1,0,
                     0,0,0,0,0,0,1),
            free=c(T, F, F, F, F, F, F,
                   F, T, F, F, F, F, F,
                   F, F, T, F, F, F, F,
                   F, F, F, T, F, F, F,
                   F, F, F, F, T, F, F,
                   F, F, F, F, F, T, F,
                   F, F, F, F, F, F, T),
            labels=c("e1", NA,   NA,   NA,   NA,   NA,   NA,
                     NA, "e2",   NA,   NA,   NA,   NA,   NA,
                     NA,   NA, "e3",   NA,   NA,   NA,   NA,
                     NA,   NA,   NA, "e4",   NA,   NA,   NA,
                     NA,   NA,   NA,   NA, "e5",   NA,   NA,
                     NA,   NA,   NA,   NA,   NA, "e6",   NA,
                     NA,   NA,   NA,   NA,   NA,   NA, "varF1"),
            byrow=TRUE,
            name="S"),
      mxMatrix("Full",  nrow=6, ncol=7,
            free=FALSE,
            values=c(1,0,0,0,0,0,0,
                     0,1,0,0,0,0,0,
                     0,0,1,0,0,0,0,
                     0,0,0,1,0,0,0,
                     0,0,0,0,1,0,0,
                     0,0,0,0,0,1,0),
            byrow=TRUE,
            dimnames=list(
            	c("x1","x2","x3","x4","x5","x6"),
            	c("x1","x2","x3","x4","x5","x6","F1")),
            name="F"),
      mxMatrix("Full", nrow=1, ncol=7,
            values=c(1,1,1,1,1,1,0),
            free=c(T,T,T,T,T,T,F),
            labels=c("meanx1","meanx2","meanx3",
                     "meanx4","meanx5","meanx6",
                     NA),
            dimnames=list(NULL, c("x1","x2","x3","x4","x5","x6","F1")),
            name="M"),
      mxRAMObjective("A","S","F","M")
      )
      
oneFactorFit<-mxRun(oneFactorModel)

omxCheckCloseEnough(oneFactorFit@output$estimate[["l2"]], 0.999, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l3"]], 0.959, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l4"]], 1.028, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l5"]], 1.008, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["l6"]], 1.021, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["varF1"]], 0.645, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e1"]], 0.350, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e2"]], 0.379, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e3"]], 0.389, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e4"]], 0.320, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e5"]], 0.370, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["e6"]], 0.346, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx1"]], 2.988, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx2"]], 3.011, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx3"]], 2.986, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx4"]], 3.053, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx5"]], 3.016, 0.01)
omxCheckCloseEnough(oneFactorFit@output$estimate[["meanx6"]], 3.010, 0.01)
