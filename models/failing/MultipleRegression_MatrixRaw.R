myRegDataRaw<-read.table("myRegData.txt",header=TRUE)

myRegDataRaw<-myRegDataRaw[,c("x","y","z")]

model<-mxModel("Multiple Regression - Matrix Specification", 
      mxData(myRegDataRaw,type="raw"),
      mxMatrix("Full", nrow=3, ncol=3,
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
            name="A"),
      mxMatrix("Symm", nrow=3, ncol=3, 
            values=c(1, 0, .5,
                     0, 1, 0,
                    .5, 0, 1),
            free=c(T, F, T,
                   F, T, F,
                   T, F, T),
            labels=c("varx", NA,         "covxz",
                     NA,     "residual",   NA,
                     "covxz",   NA,         "varz"),
            byrow=TRUE,
            name="S"),
      mxMatrix("Iden",  nrow=3, ncol=3,
            name="F"),
      mxMatrix("Full", nrow=1, ncol=3,
            values=c(0,0,0),
            free=c(T,T,T),
            labels=c("meanx","beta0","meanz"),
            name="M"),
            mxRAMObjective("A","S","F","M")
      )
      
multipleRegMatrixRaw<-mxRun(model)

multipleRegMatrixRaw@output

omxCheckCloseEnough(multipleRegMatrixRaw@output$estimate[["beta0"]], 1.6331, 0.001)
omxCheckCloseEnough(multipleRegMatrixRaw@output$estimate[["betax"]], 0.4246, 0.001)
omxCheckCloseEnough(multipleRegMatrixRaw@output$estimate[["betaz"]], 0.2260, 0.001)
omxCheckCloseEnough(multipleRegMatrixRaw@output$estimate[["residual"]], 0.646, 0.001)
omxCheckCloseEnough(multipleRegMatrixRaw@output$estimate[["varx"]], 1.116, 0.001)
omxCheckCloseEnough(multipleRegMatrixRaw@output$estimate[["varz"]], 0.836, 0.001)
omxCheckCloseEnough(multipleRegMatrixRaw@output$estimate[["cov"]], 0.289, 0.001)
omxCheckCloseEnough(multipleRegMatrixRaw@output$estimate[["meanx"]], 0.054, 0.001)
omxCheckCloseEnough(multipleRegMatrixRaw@output$estimate[["meanz"]], 4.061, 0.001)