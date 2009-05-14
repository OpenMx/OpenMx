data<-mxData(matrix(c(3.6,2.2,0.5,2.2,3.2,2,0.5,2,3.9),nrow=3),
	type="cov",
	numObs=100)

s<-mxMatrix("Symm", free=TRUE, 
	values=matrix(c(10,9,9,9,10,9,9,9,10),nrow=3),
	matrix(c("v1", "c12", "c13", 
		"c12", "v2", "c23", 
		"c13", "c23", "v3"),nrow=3),
	name="s"
	)
	
c<-mxMatrix("Full", free=FALSE,
	values=matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),nrow=3),
	name="c"
	)
	
bounds<-mxBounds(c("v1", "v2", "v3"), min=0)

constraint<-mxConstraint("s", ">", "c", name="constraint")

model<-mxModel("model", data, s, c, 
	constraint, bounds, mxMLObjective("s"))
	
run<-mxRun(model)
#Error in mxRun(model) : Non-positive-definite.