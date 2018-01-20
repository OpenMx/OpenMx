#
#   Copyright 2007-2018 The OpenMx Project
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


require(OpenMx)
myDataRaw<-read.table("data/myAutoregressiveData.txt",header=T)

model<-mxModel("Autoregressive Model, Matrix Specification, Raw Data",
      mxData(myDataRaw,type="raw"),
      mxMatrix("Full", nrow=5, ncol=5,
            values=c(0,1,0,0,0,
                     0,0,1,0,0,
                     0,0,0,1,0,
                     0,0,0,0,1,
                     0,0,0,0,0),
            free=c(F, T, F, F, F,
                   F, F, T, F, F,
                   F, F, F, T, F,
                   F, F, F, F, T,
                   F, F, F, F, F),
            labels=c(NA, "beta", NA,    NA,    NA,
                     NA, NA,    "beta", NA,    NA,
                     NA, NA,     NA,   "beta", NA,
                     NA, NA,     NA,    NA,   "beta",
                     NA, NA,     NA,    NA,    NA),
            byrow=FALSE,
            name="A"),
      mxMatrix("Symm", nrow=5, ncol=5, 
            values=c(1, 0, 0, 0, 0,
                     0, 1, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 1, 0,
                     0, 0, 0, 0, 1),
            free=c(T, F, F, F, F,
                   F, T, F, F, F,
                   F, F, T, F, F,
                   F, F, F, T, F,
                   F, F, F, F, T),
            labels=c("varx", NA,  NA,  NA,  NA,
                     NA,     "e2", NA,  NA,  NA,
                     NA,      NA, "e3", NA,  NA,
                     NA,      NA,  NA, "e4", NA,
                     NA,      NA,  NA,  NA, "e5"),
            byrow=TRUE,
            name="S"),
      mxMatrix("Iden",  nrow=5, ncol=5,
      		dimnames=list(
				c("x1","x2","x3","x4","x5"), c("x1","x2","x3","x4","x5")),
            name="F"),
      mxMatrix("Full", nrow=1, ncol=5,
            values=c(1,1,1,1,1),
            free=c(T, T, T, T, T),
            labels=c("mean1","int2","int3","int4","int5"),
            dimnames=list(
				NULL, c("x1","x2","x3","x4","x5")),
            name="M"),
      mxFitFunctionML(),mxExpectationRAM("A","S","F","M")
      )
      
autoregressiveMatrixRaw<-mxRun(model)

autoregressiveMatrixRaw$output



# Comparing to old Mx Output
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["beta"]], 0.4267, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["varx"]], 0.6658, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["e2"]], 1.1420, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["e3"]], 1.0379, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["e4"]], 0.7908, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["e5"]], 0.8176, 0.01)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["mean1"]], 3.0537, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["int2"]], 0.0824, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["int3"]], 0.0890, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["int4"]], -0.0362, 0.001)
omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["int5"]], -0.1356, 0.001)

# Comparing to Mplus values
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["beta"]], 0.427, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["varx"]], 0.665, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["e2"]], 1.142, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["e3"]], 1.038, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["e4"]], 0.791, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["e5"]], 0.818, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["mean1"]], 3.054, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["int2"]], 0.082, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["int3"]], 0.089, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["int4"]], -0.036, 0.001)
# omxCheckCloseEnough(autoregressiveMatrixRaw$output$estimate[["int5"]], -0.135, 0.001)
