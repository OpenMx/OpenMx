require(OpenMx)

data(myFADataRaw)

twoFactorRaw <- myFADataRaw[,c("x1","x2","x3","y1","y2","y3")]
manifest <- names(twoFactorRaw)
latent = c("F1", "F2")
vars <- c(manifest, latent)

twoFactorModel <- mxModel("Two Factor Model - Matrix", 
    type="RAM",
    mxData(
        observed=twoFactorRaw, 
        type="raw",
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
        name="F", dimnames=list(manifest, vars)),
    mxMatrix("Full", nrow=8, ncol=1,
        values=c(1,1,1,1,1,1,0,0),
        free=c(T,T,T,T,T,T,F,F),
        labels=c("meanx1","meanx2","meanx3",
                 "meanx4","meanx5","meanx6",
                  NA,NA),
        name="M"),
    mxRAMObjective("A","S","F","M")
)
      
twoFactorFit <- mxRun(twoFactorModel)

# Old Mx Values
omxCheckCloseEnough(twoFactorFit@output$estimate[["l2"]], 0.9723, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l3"]], 0.9313, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l5"]], 1.0498, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["l6"]], 1.0531, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["varF1"]], 0.6604, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["varF2"]], 0.4505, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["cov"]], 0.2952, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e1"]], 0.3349, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e2"]], 0.3985, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e3"]], 0.4091, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e4"]], 0.5404, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e5"]], 0.4809, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["e6"]], 0.5571, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx1"]], 2.988, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx2"]], 3.0113, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx3"]], 2.9861, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx4"]], 2.9554, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx5"]], 2.9562, 0.01)
omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx6"]], 2.9673, 0.01)

# omxCheckCloseEnough(twoFactorFit@output$estimate[["l2"]], 0.999, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["l3"]], 0.959, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["l4"]], 1.028, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["l5"]], 1.008, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["l6"]], 1.021, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["varF1"]], 0.645, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["e1"]], 0.350, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["e2"]], 0.379, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["e3"]], 0.389, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["e4"]], 0.320, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["e5"]], 0.370, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["e6"]], 0.346, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx1"]], 2.988, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx2"]], 3.011, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx3"]], 2.986, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx4"]], 3.053, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx5"]], 3.016, 0.01)
# omxCheckCloseEnough(twoFactorFit@output$estimate[["meanx6"]], 3.010, 0.01)

