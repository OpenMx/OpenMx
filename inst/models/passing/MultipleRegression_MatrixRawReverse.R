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


# SCRIPT: MultipleRegression_MatrixRawReverse.R
# Author: 
# History:  Sat 26 Sep 2009 14:57:38 BST
# (tb) updated to use data(), addded var name array, added summary statement... 
# OpenMx: http://openmx.ssri.psu.edu
##########################################
require(OpenMx)
data("myRegDataRaw", package="OpenMx")
MultipleDataRaw<-myRegData[,c("z","y","x")]
# selVars is the reverse of the selection order to test
# dimnames on the means and covariance matrices
selVars <- c("x","y","z")
multiRegModel<-mxModel("Multiple Regression Matrix Specification", 
    mxData(MultipleDataRaw,type="raw"),
    mxMatrix("Full", nrow=3, ncol=3,
        values=c(0, 0, 0,
                 1, 0, 1,
                 0, 0, 0),
        free=  c(F, F, F,
                 T, F, T,
                 F, F, F),
        labels=c(NA,     NA, NA,
                "betax", NA,"betaz",
                 NA,     NA, NA),
        byrow=TRUE,
        name = "A" ),
    mxMatrix("Symm", nrow=3, ncol=3, 
        values=c(1, 0, .5,
                 0, 1, 0,
                .5, 0, 1),
        free=  c(T, F, T,
                 F, T, F,
                 T, F, T),
        labels=c("varx",  NA,         "covxz",
                  NA,    "residual",   NA,
                 "covxz", NA,         "varz"),
        byrow=TRUE,
        name="S" ),
    mxMatrix("Iden",  nrow=3, ncol=3,
        dimnames = list(selVars,selVars),
		name="F"
	),
    mxMatrix("Full", nrow=1, ncol=3, values=0, free  =T,
        labels=c("meanx","beta0","meanz"),
        dimnames = list(NULL, selVars),
        name="M"
	),
    mxFitFunctionML(),mxExpectationRAM("A","S","F","M")
)
      
multiRegFit<-mxRun(multiRegModel)
summary(multiRegFit)

# Old Mx Output
omxCheckCloseEnough(multiRegFit$output$estimate[["beta0"]], 1.6332, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["betax"]], 0.4246, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["betaz"]], 0.2260, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["residual"]], 0.6267, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["varx"]], 1.1053, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["varz"]], 0.8275, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["covxz"]], 0.2862, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["meanx"]], 0.0542, 0.001)
omxCheckCloseEnough(multiRegFit$output$estimate[["meanz"]], 4.0611, 0.001)
