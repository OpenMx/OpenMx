#Definition Variable Test 1
#Author: Ryne Estabrook
#Date: 12 May 2009


#This script is used to test the definition variable functionality in OpenMx
#The definition variable in this example is dichotomous, and describes two different groups
#These two groups are measured on two variables, x and y
#The group with a definition value of 1 has a perfect relationship between x and y
#The group with a definition value of 0 has no relationship between x and y (beyond chance)
#The definition variable will then be used to define the covarinace between x and y


#make some data!
#two groups; in group 1, x and y are perfectly correlated
x1<-rnorm(50)
y1<-x1

#in group 0, x and y have no relationship
x2<-rnorm(50)
y2<-rnorm(50)

#put them both together, add a definition variable, and make an MxData object
x<-c(x1,x2)
y<-c(y1,y2)
def<-rep(c(1,0),each=50)
data<-mxData(data.frame(x,y,def), type="raw")

#define the model: we'll just use an S matrix and let A and F drop out
#as currently specified, this would fit a zero df model to a 2x2 covariance matrix
S<-mxMatrix("Symm", values=c(1,.5,1), free=TRUE, nrow=2, ncol=2, name="S")

#define the model, including a FIML objective function, which will optimize the matrix S
model<-mxModel("model", mxFIMLObjective("S"), data, S)

#specify a definition variable
#as.list function is required
model@objective@definitionVars<-as.list("def")

#run the model, which does not include the definition variable
#(presumably, specifying "def" as a definition variable removes it from the covariance algebra)
run<-mxRun(model)
#Error in mxRun(model) :
#REAL() can only be applied to a 'numeric', not a 'list'
#Error in args(REAL) : no function to return from, jumping to top level


#This doesn't work, but would specifying the covariance between x and y to be equal to def look
#something like this?
model[["S"]]@values[2,1]<-model@data$def
model[["S"]]@values[1,2]<-model@data$def

#or, with the new namespace ideas
model[["S"]]@values[2,1]<-model.data.def
model[["S"]]@values[1,2]<-model.data.def

#Questions:
#1. Does specifying the variable "def" in the "definitionVars" slot define a definition variable,
#    and remove it entirely from covariance algebra?
#2. Are the @ symbols the best way to access this slot?
#3. Despite the specification of type="Symm", the values matrix requires that both off-diagional positions
#    on the S matrix be changed. Is this right?