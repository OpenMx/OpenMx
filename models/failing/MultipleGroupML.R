data1<-mxData(matrix(1), type="cov", numObs=100)
data2<-mxData(matrix(2), type="cov", numObs=100)

mat1<-mxMatrix("Full",1.1,free=TRUE, nrow=1, ncol=1, name="mat1")
mat2<-mxMatrix("Full",1.9,free=TRUE, nrow=1, ncol=1, name="mat2")

obj1<-mxMLObjective("mat1")
obj2<-mxMLObjective("mat2")

model1<-mxModel("first", data1, mat1, obj1)
model2<-mxModel("second", data2, mat2, obj2)

output1<-mxRun(model1)
output2<-mxRun(model2)

output1@output
output2@output

alg<-mxAlgebra(model1.objective + model2.objective, name="alg")
obj<-mxAlgebraObjective("alg")

model<-mxModel("both", alg, obj, model1, model2)
output<-mxRun(model)
