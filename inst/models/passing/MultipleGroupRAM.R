#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
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


#This one has comments!
#Author: Ryne Estabrook
#Created: 30 Apr 2009

require(OpenMx)

#Data: 1x1 "covariance" matrices (Ok, variance matrices)
data1 <- mxData(matrix(1, dimnames = list('a', 'a')), type="cov", numObs=100)
data2 <- mxData(matrix(2, dimnames = list('a', 'a')), type="cov", numObs=100)

#S Matrices: 1 x 1 with a free parameter
S1 <- mxMatrix("Full", 1.1,free=TRUE, nrow=1, ncol=1, name="S")
S2 <- mxMatrix("Full", 2.1,free=TRUE, nrow=1, ncol=1, name="S")

#A Matrix, 1 x 1 Zero Matrix
matrixA <- mxMatrix("Zero", nrow=1, ncol=1, name="A")

#F Matrix, 1 x 1 Identity Matrix
matrixF <- mxMatrix("Iden", nrow=1, name="F")

#Lets make some expectation functions!
expectation <- mxExpectationRAM("A", "S", "F", dimnames = c('a'))

#Models
model1 <- mxModel("first", matrixA, S1, matrixF, expectation, data1, mxFitFunctionML())
model2 <- mxModel("second", matrixA, S2, matrixF, expectation, data2, mxFitFunctionML())

#Run them
output1 <- mxRun(model1, suppressWarnings = TRUE)
output2 <- mxRun(model2, suppressWarnings = TRUE)

###Starting the "Super" Model, which contains models 1 and 2
#This will use the mxFitFunctionMultigroup function

#now the objective/fit function
obj<-mxFitFunctionMultigroup(c('first', 'second'))

#make a model
model<-mxModel("both", obj, model1, model2)

#run the "super" model
output<-mxRun(model, suppressWarnings = TRUE)

###Check Results
#Model 1: This should have a value of 1
coef(output1)

#Model 2: This should have a value of 2
coef(output2)

#"Super" Model: This should have values of 1 and 2, in order.
coef(output)

omxCheckCloseEnough(coef(output1), .99 * c(1), 0.001)
omxCheckCloseEnough(coef(output2), .99 * c(2), 0.001)
omxCheckCloseEnough(coef(output), .99 * c(1, 2), 0.001)

