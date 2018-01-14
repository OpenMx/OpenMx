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


#Amended version of MultipleGroupRAMconstraint.R, with constraint on free parameters
#Author: Ryne Estabrook
#Created: 30 Apr 2009

#Goal: Constrain the single parameter in each group to be equal

require(OpenMx)

#Data: 1x1 "covariance" matrices (Ok, variance matrices)
data1 <- mxData(matrix(1, dimnames = list('a', 'a')), type="cov", numObs=100)
data2 <- mxData(matrix(2, dimnames = list('a', 'a')), type="cov", numObs=100)

#S Matrices: 1 x 1 with a free parameter (must have the same value in multiple group estimation)
S1 <- mxMatrix("Full", 1.5, free=TRUE, nrow=1, ncol=1, labels="parameter", name="S")
S2 <- mxMatrix("Full", 1.5, free=TRUE, nrow=1, ncol=1, labels="parameter", name="S")

#A Matrix, 1 x 1 Zero Matrix
matrixA <- mxMatrix("Zero", nrow=1, ncol=1, name="A")

#F Matrix, 1 x 1 Identity Matrix
matrixF <- mxMatrix("Iden", nrow=1, name="F", dimnames = list('a', 'a'))

#Lets make some objective functions!
objective <- mxExpectationRAM("A", "S", "F")

#Models
model1<-mxModel("first", matrixA, S1, matrixF, objective, data1, mxFitFunctionML())
model2<-mxModel("second", matrixA, S2, matrixF, objective, data2, mxFitFunctionML())

#Run them
output1<-mxRun(model1, suppressWarnings=TRUE)
output2<-mxRun(model2, suppressWarnings=TRUE)

###Starting the "Super" Model, which contains models 1 and 2
#This will use the mxFitFunctionAlgebra function
#we first need an algebra, which describes how obj1 and obj2 go together (sum)
alg<-mxAlgebra(first.objective + second.objective, name="alg")

#now the objective function
obj <- mxFitFunctionAlgebra("alg")

#make a model
model <- mxModel("both", alg, obj, model1, model2)

#run the "super" model
output<-mxRun(model, suppressWarnings=TRUE)

###Check Results
#Model 1: This should have a value of 1
print(output1$output$estimate)

#Model 2: This should have a value of 2
print(output2$output$estimate)

#"Super" Model: This should have a value of 1.5
print(output$output$estimate)

omxCheckCloseEnough(output1$output$estimate, .99 * c(1), 0.001)
omxCheckCloseEnough(output2$output$estimate, .99 * c(2), 0.001)
omxCheckCloseEnough(output$output$estimate, .99 * c(1.5), 0.001)
