#Amended version of MultipleGroupRAMconstraint.R, with constraint on free parameters
#Author: Ryne Estabrook
#Created: 30 Apr 2009

#Goal: Constrain the single parameter in each group to be equal

#Data: 1x1 "covariance" matrices (Ok, variance matrices)
data1<-mxData(matrix(1), type="cov", numObs=100, name="data1")
data2<-mxData(matrix(2), type="cov", numObs=100, name="data2")

#S Matrices: 1 x 1 with a free parameter (must have the same value in multiple group estimation)
mat1<-mxMatrix("Full",1.5,free=TRUE, nrow=1, ncol=1, labels="m1", name="mat1")
mat2<-mxMatrix("Full",1.5,free=TRUE, nrow=1, ncol=1, labels="m1", name="mat2")

#A Matrix, 1 x 1 Zero Matrix
#If the same A matrix is used in both groups, mxRun throws and error.
a1<-mxMatrix("Zero", nrow=1, ncol=1, name="A1")
a2<-mxMatrix("Zero", nrow=1, ncol=1, name="A2")

#F Matrix, 1 x 1 Identity Matrix
#Same error as the A matrices
f1<-mxMatrix("Iden", nrow=1, name="F1")
f2<-mxMatrix("Iden", nrow=1, name="F2")

#Lets make some objective functions!
obj1<-mxRAMObjective("A1","mat1","F1",name="obj1")
obj2<-mxRAMObjective("A2","mat2","F2", name="obj2")

#Models
model1<-mxModel("first",data1, mat1, a1, f1, obj1)
model2<-mxModel("second",data2, mat2, a2, f2, obj2)

#Run them
output1<-mxRun(model1)
output2<-mxRun(model2)

###Starting the "Super" Model, which contains models 1 and 2
#This will use the mxAlgebraObjective function
#we first need an algebra, which describes how obj1 and obj2 go together (sum)
alg<-mxAlgebra(obj1+obj2, name="alg")

#now the objective function
obj<-mxAlgebraObjective("alg", name="obj")

#make a model
model<-mxModel("both", alg, obj, model1, model2)

#run the "super" model
output<-mxRun(model)

###Check Results
#Model 1: This should have a value of 1
output1@output$estimate

#Model 2: This should have a value of 2
output2@output$estimate

#"Super" Model: This should have a value of 1.5
output@output$estimate


#Notes:
#  Each model requires its own unique set of matrices:
#       if one tries to share the A and F matrices, mxRun throws an error
#  mxFIMLObjective is not as precise as it needs to be

#Related scripts:
#  MultipleGroupML.R is failing, due to mxMLObjective problems.
#  MultipleGroupRAM.R works, withing precision of mxRAMObjective